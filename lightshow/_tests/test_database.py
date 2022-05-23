from pathlib import Path
from tempfile import TemporaryDirectory

from lightshow.database import Database
from lightshow.parameters.feff import FEFFParameters


class TestDatabase:
    @staticmethod
    def test_from_materials_project(mp_Structure_mp390, mp_Structure_mvc11115):
        dat = Database.from_materials_project(["mp-390", "mvc-11115"])
        assert set(dat.structures.keys()) == {"mp-390", "mvc-11115"}
        assert dat.structures["mp-390"] == mp_Structure_mp390
        assert dat.structures["mvc-11115"] == mp_Structure_mvc11115

    @staticmethod
    def test_mpids_list(mp_Structure_mp390, mp_Structure_mvc11115):
        with TemporaryDirectory() as tempDir:

            # Create a temporary directory structure and save the structures
            # as CONTCAR files
            root = Path(tempDir) / Path("iama") / Path("test")
            path1 = root / Path("d1")
            path1.mkdir(exist_ok=True, parents=True)
            path1 = path1 / Path("CONTCAR")
            mp_Structure_mp390.to(fmt="POSCAR", filename=path1)
            path2 = root / Path("d2")
            path2.mkdir(exist_ok=True, parents=True)
            path2 = path2 / Path("CONTCAR")
            mp_Structure_mvc11115.to(fmt="POSCAR", filename=path2)

            dat = Database.from_files(root, filename="CONTCAR")
            assert str(path1.parent) in dat.structures.keys()
            assert str(path2.parent) in dat.structures.keys()

    @staticmethod
    def test_write(TiO10):

        # Write a few test files to disk somewhere, then read them in
        with TemporaryDirectory() as tempDir:
            root = Path(tempDir) / Path("iama") / Path("test")
            for key, struct in TiO10.items():
                path = root / Path(f"{key}")
                path.mkdir(exist_ok=True, parents=True)
                path = path / Path("POSCAR")
                struct.to(fmt="POSCAR", filename=path)

            # Load it all in
            dat = Database.from_files(root, filename="POSCAR")
            root = Path(tempDir) / Path("iama") / Path("destination")

            # Write it (note these are Ti-O compounds)
            feff_parameters = FEFFParameters(
                cards={
                    "S02": "0",
                    "COREHOLE": "RPA",
                    "CONTROL": "1 1 1 1 1 1",
                    "XANES": "4 0.04 0.1",
                    "SCF": "7.0 0 100 0.2 3",
                    "FMS": "9.0 0",
                    "EXCHANGE": "0 0.0 0.0 2",
                    "RPATH": "-1",
                },
                edge="K",
                radius=9.0,
                spectrum="XANES",
                name="FEFF",
            )
            dat.write(root, absorbing_atom="Ti", options=[feff_parameters])

from pathlib import Path
from tempfile import TemporaryDirectory

from lightshow.database import Database


class TestDatabase:
    @staticmethod
    def test_from_materials_project(Structure_mp390, Structure_mvc11115):
        dat = Database.from_materials_project(["mp-390", "mvc-11115"])
        assert set(dat.structures.keys()) == {"mp-390", "mvc-11115"}
        assert dat.structures["mp-390"] == Structure_mp390
        assert dat.structures["mvc-11115"] == Structure_mvc11115

    @staticmethod
    def test_mpids_list(Structure_mp390, Structure_mvc11115):
        with TemporaryDirectory() as tempDir:

            # Create a temporary directory structure and save the structures
            # as CONTCAR files
            root = Path(tempDir) / Path("iama") / Path("test")
            path1 = root / Path("d1")
            path1.mkdir(exist_ok=True, parents=True)
            path1 = path1 / Path("CONTCAR")
            Structure_mp390.to(fmt="POSCAR", filename=path1)
            path2 = root / Path("d2")
            path2.mkdir(exist_ok=True, parents=True)
            path2 = path2 / Path("CONTCAR")
            Structure_mvc11115.to(fmt="POSCAR", filename=path2)

            dat = Database.from_files(root, filename="CONTCAR")
            assert str(path1.parent) in dat.structures.keys()
            assert str(path2.parent) in dat.structures.keys()

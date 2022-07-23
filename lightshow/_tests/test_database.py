import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from pymatgen.core.structure import IStructure
from pymatgen.ext.matproj import MPRestError
from warnings import warn

from lightshow.database import Database
from lightshow.parameters.feff import FEFFParameters
from lightshow.parameters.vasp import VASPParameters, Incar


def _validate_single_VASP_FEFF_calculation(
    feff_input_path, vasp_input_path, r=4.0
):
    """Executes a validation scheme for testing if FEFF and VASP
    calculations produced input files in which the indexes correspond.

    Parameters
    ----------
    directory : str
        The directory of the input files to check. Usually this is an mpid.
        The directory must contain a "VASP" subdirectory and a
        "FEFF-XANES" subdirectory.
    """

    # Load FEFF
    with open(feff_input_path, "r") as f:
        feff_lines = f.readlines()
    where_atoms = [
        ii for ii, xx in enumerate(feff_lines) if xx.startswith("ATOMS")
    ]
    assert len(where_atoms) == 1
    feff_lines = feff_lines[where_atoms[0] + 3 : -1]
    feff_lines = [xx.split() for xx in feff_lines]

    # Exclude self
    feff_distances = (
        np.array([xx[5] for xx in feff_lines[1:]]).astype(float).round(6)
    )

    # Load VASP
    with open(vasp_input_path, "r") as f:
        vasp_lines = f.readlines()

    # VASP POSCAR files are easy, only need data after line 8
    vasp_lines = [xx.strip().split() for xx in vasp_lines[8:-1]]
    vasp_coordinates = np.array([xx[:3] for xx in vasp_lines]).astype(float)
    vasp_atoms = [xx[3] for xx in vasp_lines]

    # Load in the pymatgen structure representation of the VASP POSCAR
    vasp_structure = IStructure.from_file(vasp_input_path)

    # Check site ordering is the same
    assert [site.specie.symbol for site in vasp_structure] == vasp_atoms

    # Sanity check that the coordinates are the same
    vasp_coordinates_from_struct = np.array(
        [site.frac_coords for site in vasp_structure]
    )
    np.testing.assert_array_almost_equal(
        vasp_coordinates_from_struct, vasp_coordinates
    )

    # The FEFF atom ordering and coordinates are not the same. However,
    # the nearest neighbor distance should be. The FEFF files are output
    # in order of increasing distances from the first atom.
    neigh = vasp_structure.get_neighbors(vasp_structure[0], r=r)
    distances = np.array([xx.nn_distance for xx in neigh])
    distances.sort()
    distances = distances.round(6)
    min_dist = min(len(distances), len(feff_distances))

    return distances[:min_dist], feff_distances[:min_dist]


def _validate_VASP_FEFF_calculations_match(root, decimals=4, r=6.0):
    """Validates all XAS input files against FEFF input files in a
    directory. See _validate_single_VASP_FEFF_calculation for more details.

    Parameters
    ----------
    decimals : int, optional
        The number of decimals to validate the structures up to. Default is
        5. Recommended to not change this, since it corresponds somewhat
        to the number of decimals saved to disk.
    r : float, optional
        The radius up to which to construct/consider nearest neighbors for
        checking the bond distances.
    verbose : bool, optional
        If True, outputs a progress bar
    """

    directories = [xx for xx in list(Path(root).iterdir()) if xx.is_dir()]

    for path in directories:

        # Try the XANES directory by default
        p1 = path / Path("FEFF")
        p2 = path / Path("VASP")

        for subdir in p1.iterdir():
            feff_input_path = Path(subdir) / Path("feff.inp")
            vasp_input_path = Path(p2) / Path(subdir.name) / Path("POSCAR")

            d1, d2 = _validate_single_VASP_FEFF_calculation(
                feff_input_path, vasp_input_path, r=r
            )

            np.testing.assert_array_almost_equal(d1, d2, decimal=decimals)


class TestDatabase:
    @staticmethod
    def test_from_materials_project(mp_Structure_mp390, mp_Structure_mvc11115):
        try:
            dat = Database.from_materials_project(
                ["mp-390", "mvc-11115"], suppress_MPRestError=False
            )
        except MPRestError:
            warn(
                "MPRestError during pulling data- using locally stored "
                "Database for the test"
            )
            path = Path.cwd() / Path("lightshow") / Path("_tests")
            path = path / Path("database_obj_test.json")
            with open(path, "r") as infile:
                loaded_json = json.load(infile)
            dat = Database.from_dict(loaded_json)
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
    def test_example_database(TiO10, dummy_potcar_file_directory):

        # Write a few test files to disk somewhere, then read them in
        with TemporaryDirectory() as tempDir:
            root = Path(tempDir) / Path("iama") / Path("test")
            for key, struct in TiO10.items():
                path = root / Path(f"{key}")
                path.mkdir(exist_ok=True, parents=True)
                path = path / Path("POSCAR")
                struct.to(fmt="POSCAR", filename=path)
                print("wrote to ", path)

            # Load it all in
            dat = Database.from_files(root, filename="POSCAR")

            target = Path(tempDir) / Path("iama") / Path("destination")
            target.mkdir(exist_ok=True, parents=True)

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
                radius=6.0,
                spectrum="XANES",
                name="FEFF",
            )
            vasp_params_corehole = VASPParameters(
                incar=Incar.from_default(neutral=False),
                potcar_directory=dummy_potcar_file_directory,
                force_spin_unpolarized=False,
            )
            dat.write(
                target,
                absorbing_atom="Ti",
                options=[feff_parameters, vasp_params_corehole],
            )

            _validate_VASP_FEFF_calculations_match(root)

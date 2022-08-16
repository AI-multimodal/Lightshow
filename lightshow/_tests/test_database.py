from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from pymatgen.core.structure import IStructure
from pymatgen.ext.matproj import MPRestError
from warnings import warn

from lightshow import (
    Database,
    FEFFParameters,
    VASPParameters,
    OCEANParameters,
    XSpectraParameters,
    EXCITINGParameters,
)
from lightshow.defaults import VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL


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
    def test_database_from_disk(database_from_file, test_structure_names):
        dat = database_from_file
        dat.cleanup_paths()
        dat.initialize_supercells(9.0)
        dat.initialize_inequivalent_sites()
        assert set(dat.structures.keys()) == set(test_structure_names)
        assert set(dat.metadata.keys()) == set(test_structure_names)
        assert set(dat.supercells.keys()) == set(test_structure_names)

    @staticmethod
    def test_from_materials_project(test_structure_names):
        try:
            dat = Database.from_materials_project(test_structure_names)
            dat.initialize_supercells(9.0)
            assert set(dat.structures.keys()) == set(test_structure_names)
            assert set(dat.supercells.keys()) == set(test_structure_names)
        except MPRestError:
            warn("MPRestError during pulling data")

    @staticmethod
    def test_write(dummy_potcar_file_directory, database_from_file):

        # Write a few test files to disk somewhere, then read them in
        with TemporaryDirectory() as tempDir:

            # Load it all in
            dat = database_from_file
            dat.initialize_supercells(9.0)
            dat.initialize_inequivalent_sites()

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
                incar=VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL,
                potcar_directory=dummy_potcar_file_directory,
                force_spin_unpolarized=False,
            )
            ocean_params = OCEANParameters(edge="K")
            exciting_params = EXCITINGParameters(edge="K")
            xspectra_params = XSpectraParameters(edge="K")

            dat.write(
                target,
                absorbing_atoms=["Ti", "O"],
                options=[
                    feff_parameters,
                    vasp_params_corehole,
                    ocean_params,
                    exciting_params,
                    xspectra_params,
                ],
            )

            _validate_VASP_FEFF_calculations_match(target)

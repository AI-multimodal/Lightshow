from pathlib import Path
import pytest
import sys

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

# Helper testing files
sys.path.append(str(Path(__file__).parent.resolve() / Path("helpers")))
from geometry import consistency_check  # noqa


def test_database_from_disk(database_from_file, test_structure_names):
    dat = database_from_file
    dat.initialize_supercells(9.0)
    dat.initialize_inequivalent_sites()
    assert set(dat.structures.keys()) == set(test_structure_names)
    assert set(dat.metadata.keys()) == set(test_structure_names)
    assert set(dat.supercells.keys()) == set(test_structure_names)


def test_from_materials_project(test_structure_names):
    try:
        dat = Database.from_materials_project(test_structure_names)
        dat.initialize_supercells(9.0)
        assert set(dat.structures.keys()) == set(test_structure_names)
        assert set(dat.supercells.keys()) == set(test_structure_names)
    except MPRestError:
        warn("MPRestError during pulling data")


@pytest.mark.parametrize(
    "mpid",
    [
        "mp-390",
        "mvc-11115",
        "mp-1215",
        "mp-1840",
        "mp-2657",
        "mp-2664",
        "mp-430",
        "mp-458",
        "mp-10734",
    ],
)
def test_write(
    mpid,
    dummy_potcar_file_directory,
    dummy_psp_file_directory,
    dummy_chpsp_file_directory,
    database_from_file,
    tmp_path,
):

    # Load it all in
    dat = database_from_file
    dat.initialize_supercells(9.0)
    dat.initialize_inequivalent_sites()
    dat._supercells = {mpid: dat.supercells[mpid]}
    dat._structures = {mpid: dat.structures[mpid]}
    dat._metadata = {mpid: dat.metadata[mpid]}

    target = Path(tmp_path) / Path("iama") / Path("destination")
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
    xspectra_params = XSpectraParameters(
        psp_directory=dummy_psp_file_directory,
        psp_cutoff_table="mock_cutoff_table.json",
        chpsp_directory=dummy_chpsp_file_directory,
        edge="K",
    )

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

    # Assert geometires
    consistency_check(target / Path(mpid), rounding=3)

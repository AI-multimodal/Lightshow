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
from lightshow._tests.conftest import get_mpids_for_stress_test
from lightshow.common.kpoints import GenericEstimatorKpoints
from lightshow.common.nbands import UnitCellVolumeEstimate

# Helper testing files
sys.path.append(str(Path(__file__).parent.resolve() / Path("helpers")))
from geometry import consistency_check  # noqa  # type: ignore


def test_from_materials_project(test_from_materials_project_structure_names):
    try:
        dat = Database.from_materials_project(
            material_ids=test_from_materials_project_structure_names
        )
        dat.initialize_supercells(9.0)
        assert set(dat.structures.keys()) == set(
            test_from_materials_project_structure_names
        )
        assert set(dat.supercells.keys()) == set(
            test_from_materials_project_structure_names
        )
    except MPRestError:
        warn("MPRestError during pulling data")


@pytest.mark.parametrize("mpid", get_mpids_for_stress_test())
def test_geometry_stress(
    mpid,
    dummy_potcar_file_directory,
    dummy_psp_file_directory,
    dummy_chpsp_file_directory,
    database_for_stress_test,
    tmp_path,
):
    # Load it all in
    dat = database_for_stress_test
    dat._structures = {mpid: dat.structures[mpid]}
    dat.initialize_supercells(9.0)
    dat.initialize_inequivalent_sites()
    dat._supercells = {mpid: dat.supercells[mpid]}
    dat._metadata = {mpid: dat.metadata[mpid]}

    target = Path(tmp_path) / Path("iama") / Path("destination")
    target.mkdir(exist_ok=True, parents=True)

    R = 10.0

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
        radius=R - 1.0,
        spectrum="XANES",
        name="FEFF",
    )
    vasp_params_corehole = VASPParameters(
        incar=VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL,
        potcar_directory=None,
        force_spin_unpolarized=False,
        kpoints=GenericEstimatorKpoints(cutoff=33),
        nbands=UnitCellVolumeEstimate(e_range=40),
    )
    ocean_params = OCEANParameters(
        edge="K",
        kpoints=GenericEstimatorKpoints(cutoff=33),
        nbands=UnitCellVolumeEstimate(e_range=40),
    )
    exciting_params = EXCITINGParameters(
        edge="K",
        kpoints=GenericEstimatorKpoints(cutoff=33),
        nbands=UnitCellVolumeEstimate(e_range=40),
    )
    xspectra_params = XSpectraParameters(
        psp_directory=None,
        psp_cutoff_table="mock_cutoff_table.json",
        chpsp_directory=None,
        edge="K",
        kpoints=GenericEstimatorKpoints(cutoff=33),
    )

    try:
        dat.write(
            target,
            absorbing_atoms="all",
            options=[
                feff_parameters,
                vasp_params_corehole,
                ocean_params,
                exciting_params,
                xspectra_params,
            ],
            write_unit_cells=True,
            pbar=False,
        )

    except KeyError:
        # This happens when an element is present in the material that is
        # simply not compatible with some of the Lightshow features
        return

    # Assert geometires
    consistency_check(target / Path(mpid), R)

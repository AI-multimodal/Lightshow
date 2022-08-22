from pathlib import Path
import pytest
import sys

from lightshow import VASPParameters
from lightshow.defaults import VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL

# Helper testing files
sys.path.append(str(Path(__file__).parent.resolve() / Path("helpers")))
from geometry import consistency_check  # noqa


@pytest.mark.parametrize(
    "mpid",
    ["mp-390", "mvc-11115"],
)
def test_write(mpid, dummy_potcar_file_directory, database_from_file, tmp_path):

    # Load it all in
    dat = database_from_file
    dat.initialize_supercells(9.0)
    dat.initialize_inequivalent_sites()
    dat._supercells = {mpid: dat.supercells[mpid]}
    dat._structures = {mpid: dat.structures[mpid]}
    dat._metadata = {mpid: dat.metadata[mpid]}

    target = Path(tmp_path) / Path("iama") / Path("destination")
    target.mkdir(exist_ok=True, parents=True)

    vasp_params_corehole = VASPParameters(
        incar=VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL,
        potcar_directory="this/dir/does/not/exist",
        force_spin_unpolarized=False,
    )

    dat.write(
        target,
        absorbing_atoms=["Ti", "O"],
        options=[vasp_params_corehole],
        pbar=False,
    )

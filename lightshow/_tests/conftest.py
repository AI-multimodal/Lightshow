from copy import deepcopy
from pathlib import Path
import pytest

from pymatgen.core.structure import Structure

from lightshow import Database

STRUCTURE_FILES_PATH = (
    Path.cwd() / Path("lightshow") / Path("_tests") / Path("structure_files")
)


@pytest.fixture
def test_structure_names():
    return [
        "mp-390",
        "mvc-11115",
        "mp-1215",
        "mp-1840",
        "mp-2657",
        "mp-2664",
        "mp-430",
        "mp-458",
        "mp-10734",
    ]


@pytest.fixture
def database_from_file():
    dat = Database.from_files(
        STRUCTURE_FILES_PATH, filename="POSCAR", cleanup_paths=True
    )
    return deepcopy(dat)


@pytest.fixture
def dummy_potcar_file_directory():
    return str(Path(__file__).parent.resolve() / Path("dummy_potcar_files"))


@pytest.fixture
def dummy_psp_file_directory():
    return str(Path(__file__).parent.resolve() / Path("dummy_psp_files"))


@pytest.fixture
def dummy_chpsp_file_directory():
    return str(Path(__file__).parent.resolve() / Path("dummy_chpsp_files"))


@pytest.fixture
def mp_Structure_mp390():
    return Structure.from_file(
        STRUCTURE_FILES_PATH / Path("mp-390") / Path("POSCAR")
    )

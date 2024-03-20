import random
from copy import deepcopy
from functools import cache
from pathlib import Path

import pytest
from pymatgen.core.structure import Structure

from lightshow import Database

STRUCTURE_FILES_PATH = Path(__file__).parent / Path("structure_files")

SPECTRA_FILES_PATH = Path.cwd() / Path("notebooks") / Path("spectra_files")

SAMPLE_FILES_PATH = (
    Path.cwd() / Path("lightshow") / Path("_tests") / Path("sample_files")
)


@pytest.fixture
def test_structure_names():
    return [
        "mp-390",
        # "mvc-11115",  # No longer in v2 API
        "mp-1215",
        "mp-1840",
        "mp-2657",
        "mp-2664",
        "mp-430",
        "mp-458",
        "mp-10734",
    ]


@pytest.fixture
def test_from_materials_project_structure_names():
    return [
        "mp-390",
        "mp-1215",
        "mp-1840",
        "mp-2657",
        "mp-2664",
        "mp-430",
        "mp-458",
        "mp-10734",
        "mp-1203",
    ]


@pytest.fixture
def database_from_file():
    dat = Database.from_files(STRUCTURE_FILES_PATH, filename="POSCAR")
    return deepcopy(dat)


@cache
def _database_for_stress_test():
    # print("✅✅✅✅ Pulled some data from the Materials Project ✅✅✅✅")
    # print("✅✅✅✅ hopefully only once!!!                      ✅✅✅✅")
    # trouble = ["mp-1104067", "mp-1103465", "mp-1196561", "mp-1191564",
    # "mp-1208324", "mp-10734"]
    # return Database.from_materials_project(material_ids=trouble)
    # return Database.from_materials_project(material_ids=["mp-980945"])
    db = Database.from_materials_project(chemsys=["Ti-*", "Mn-O-*"])
    # db = Database.from_materials_project(material_ids=["mp-390"])
    keys = sorted(list(db._structures.keys()))
    random.seed(123)
    keys = random.sample(keys, 100)
    db._structures = {key: db._structures[key] for key in keys}
    db._metadata = {key: db._metadata[key] for key in keys}
    return db


@pytest.fixture
def database_for_stress_test():
    return deepcopy(_database_for_stress_test())


def get_mpids_for_stress_test():
    db = _database_for_stress_test()
    return list(db._structures.keys())


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


@pytest.fixture
def spectra_file_directory():
    return SPECTRA_FILES_PATH


@pytest.fixture
def sample_file_directory():
    return SAMPLE_FILES_PATH

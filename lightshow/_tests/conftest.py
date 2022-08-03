import json
from os import environ
from pathlib import Path

import pytest

from pymatgen.core.structure import Structure


API_KEY = environ.get("PMG_API_KEY", None)


@pytest.fixture
def mp_Structure_mp390():
    path = Path.cwd() / Path("lightshow") / Path("_tests")
    path = path / Path("mp-390.json")
    with open(path, "r") as infile:
        loaded_json = json.load(infile)
    return Structure.from_dict(loaded_json)


@pytest.fixture
def mp_Structure_mvc11115():
    path = Path.cwd() / Path("lightshow") / Path("_tests")
    path = path / Path("mvc-11115.json")
    with open(path, "r") as infile:
        loaded_json = json.load(infile)
    return Structure.from_dict(loaded_json)


@pytest.fixture
def TiO10():
    path = Path(__file__).parent.resolve() / Path("TiO10.json")
    with open(path, "r") as f:
        d = json.load(f)
    for key in d.keys():
        d[key] = Structure.from_dict(d[key])
    return d


@pytest.fixture
def dummy_potcar_file_directory():
    return str(Path(__file__).parent.resolve() / Path("dummy_potcar_files"))

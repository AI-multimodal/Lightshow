import json
from os import environ
from pathlib import Path

import pytest

from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester


API_KEY = environ.get("PMG_API_KEY", None)


@pytest.fixture
def mp_Structure_mp390():
    with MPRester(API_KEY) as mpr:
        metadata = mpr.get_doc("mp-390")
    return Structure.from_dict(metadata["structure"])


@pytest.fixture
def mp_Structure_mvc11115():
    with MPRester(API_KEY) as mpr:
        metadata = mpr.get_doc("mvc-11115")
    return Structure.from_dict(metadata["structure"])


@pytest.fixture
def TiO10():
    path = Path(__file__).parent.resolve() / Path("TiO10.json")
    with open(path, "r") as f:
        d = json.load(f)
    for key in d.keys():
        d[key] = Structure.from_dict(d[key])
    return d

from os import environ

import pytest

from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester


API_KEY = environ.get("PMG_API_KEY", None)


@pytest.fixture
def Structure_mp390():
    with MPRester(API_KEY) as mpr:
        metadata = mpr.get_doc("mp-390")
    return Structure.from_dict(metadata["structure"])


@pytest.fixture
def Structure_mvc11115():
    with MPRester(API_KEY) as mpr:
        metadata = mpr.get_doc("mvc-11115")
    return Structure.from_dict(metadata["structure"])

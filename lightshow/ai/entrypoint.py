import pathlib
import sys

from pymatgen.core import Structure as PymatgenStructure

from .models import XASModel

PARENT_DIRECTORY = pathlib.Path(__file__).parent.resolve()
XASBLOCKS_PATH = PARENT_DIRECTORY / "models" / "xasblock" / "v1.1.1"
AVAILABLE_COMBINATIONS = [f.stem for f in XASBLOCKS_PATH.glob("*.ckpt")]


def entrypoint():
    # TODO: proper CLI
    spectroscopy_type = sys.argv[1]
    el = sys.argv[2]
    path = sys.argv[3]
    struct = PymatgenStructure.from_file(path)
    site_idxs = [
        ii for ii, site in enumerate(struct.sites) if site.specie.symbol == el
    ]
    if len(site_idxs) == 0:
        raise ValueError(f"element {el} not found in provided structure")
    spec = XASModel(element=el, spectroscopy_type=spectroscopy_type).predict(
        struct
    )
    result = {ii: spec[ii] for ii in site_idxs}
    print(result)

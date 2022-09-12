from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

from .database import Database  # noqa
from .parameters.feff import FEFFParameters  # noqa
from .parameters.vasp import VASPParameters  # noqa
from .parameters.exciting import EXCITINGParameters  # noqa
from .parameters.ocean import OCEANParameters  # noqa
from .parameters.xspectra import XSpectraParameters  # noqa

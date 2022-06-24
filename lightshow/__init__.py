from os import environ

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def _get_API_key_from_environ():
    """Checks for an environment variable PMG_API_KEY. If does not exist,
    returns None.

    Returns
    -------
    str
    """

    return environ.get("PMG_API_KEY", None)


from lightshow.database import Database  # noqa
from lightshow.parameters.feff import FEFFParameters  # noqa
from lightshow.parameters.vasp import VASPParameters  # noqa
from lightshow.parameters.exciting import EXCITINGParameters  # noqa
from lightshow.parameters.ocean import OCEANParameters  # noqa
from lightshow.parameters.xspectra import XSpectraParameters  # noqa

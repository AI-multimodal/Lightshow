from os import environ

from lightshow._version import __version__
from lightshow.database import Database
from lightshow.parameters.exciting import EXCITINGParameters
from lightshow.parameters.feff import FEFFParameters
from lightshow.parameters.ocean import OCEANParameters
from lightshow.parameters.vasp import VASPParameters
from lightshow.parameters.xspectra import XSpectraParameters


def _get_API_key_from_environ():
    """Checks for an environment variable MP_API_KEY. If does not exist,
    returns None. Note that this is now for the Materials Project v2 api.

    Returns
    -------
    str
    """

    return environ.get("MP_API_KEY", None)


def _get_POTCAR_DIRECTORY_from_environ():
    """Checks for an environment variable VASP_POTCAR_DIRECTORY. If does not
    exist, returns None.

    Returns
    -------
    str
    """

    return environ.get("VASP_POTCAR_DIRECTORY", None)


def _get_CHPSP_DIRECTORY_from_environ():
    """Checks for an environment variable XS_CHPSP_DIRECTORY. If does not
    exist, returns None.

    Returns
    -------
    str
    """

    return environ.get("XS_CHPSP_DIRECTORY", None)


def _get_PSP_DIRECTORY_from_environ():
    """Checks for an environment variable XS_PSP_DIRECTORY. If does not
    exist, returns None.

    Returns
    -------
    str
    """

    return environ.get("XS_PSP_DIRECTORY", None)


def _get_SPECIES_DIRECTORY_from_environ():
    """Checks for an environment varialble SPECIES_DIRECRORY. If does not
    exist, returns None.

    Returns
    -------
    str
    """

    return environ.get("SPECIES_DIRECTORY", None)


__all__ = [
    "__version__",
    "Database",
    "EXCITINGParameters",
    "FEFFParameters",
    "OCEANParameters",
    "VASPParameters",
    "XSpectraParameters",
]

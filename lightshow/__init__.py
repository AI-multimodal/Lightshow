from os import environ


# DO NOT CHANGE BELOW ---------------------------------------------------------
# This is replaced at build time automatically during deployment and
# installation. Replacing anything will mess that up and crash the entire
# build.
# __version__ = ...  # semantic-version-placeholder
# DO NOT CHANGE ABOVE ---------------------------------------------------------

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


from .database import Database  # noqa
from .parameters.feff import FEFFParameters  # noqa
from .parameters.vasp import VASPParameters  # noqa
from .parameters.exciting import EXCITINGParameters  # noqa
from .parameters.ocean import OCEANParameters  # noqa
from .parameters.xspectra import XSpectraParameters  # noqa

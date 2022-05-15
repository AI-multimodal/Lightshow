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

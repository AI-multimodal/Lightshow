from copy import copy
from pathlib import Path
from warnings import warn

from monty.json import MSONable
from pymatgen.io.feff.sets import MPXANESSet, MPEXAFSSet, FEFFDictSet

from lightshow.parameters._base import _BaseParameters


FEFF_DEFAULT_CARDS = {
    "S02": "0",
    "COREHOLE": "RPA",
    "CONTROL": "1 1 1 1 1 1",
    "XANES": "4 0.04 0.1",
    "SCF": "7.0 0 100 0.2 3",
    "FMS": "9.0 0",
    "EXCHANGE": "0 0.0 0.0 2",
    "RPATH": "-1",
}


class FEFFParameters(MSONable, _BaseParameters):
    """A one-stop-shop for all the different ways to modify input parameters
    for a FEFF calculation. This class is a lightweight wrapper for the
    ``FEFFDictSet``, containing some extra information and methods used for
    writing the appropriate input files. Like most Pymatgen objects, the
    class is serializable via e.g. ``feff_parameters.as_dict()``.

    Parameters
    ----------
    cards : dict
        A dictionary of of the cards to be used in the FEFF calculations.
        For example, for XANES, one might wish to use something like

        .. code-block:: python

            cards = {
                "S02": "0",
                "COREHOLE": "RPA",
                "CONTROL": "1 1 1 1 1 1",
                "XANES": "4 0.04 0.1",
                "SCF": "7.0 0 100 0.2 3",
                "FMS": "9.0 0",
                "EXCHANGE": "0 0.0 0.0 2",
                "RPATH": "-1"
            }

        And for EXAFS,

        .. code-block:: python

            cards = {
                "S02": "0",
                "COREHOLE": "RPA",
                "CONTROL": "1 1 1 1 1 1",
                "EXAFS": "20",
                "RPATH": "10",
                "SCF": "7.0 0 100 0.2 3",
                "EXCHANGE": "0 0.0 0.0 0"
            }
    edge : str
        The XAS edge of the calculation.
    radius : float
        FEFF uses clusters for its calculations. The ``radius`` parameter
        determines how large to make the cluster. It is calculated from the
        absorbing atom center in units of Angstroms.
    spectrum : str
        The type of spectroscopy to be run. This is used to set the default
        cards through Pymatgen. These defaults can be overridden by setting
        the ``cards`` argument when instantiating the class. The type of
        calculation to run. Should likely be either ``"XANES"`` or ``"EXAFS"``.
    **feff_dict_set_kwargs
        Keyword arguments to pass directly to the FEFFDictSet object before
        writing the input files. For example, ``nkpts`` can be used to specify
        the number of k-points used in the Brillouin zone. It is only used if
        FEFF is run in reciprocal space mode.
    """

    def __init__(
        self,
        cards=FEFF_DEFAULT_CARDS,
        edge="K",
        radius=9.0,
        spectrum="XANES",
        name="FEFF",
        **feff_dict_set_kwargs,
    ):
        # Try to see if "edge" is in the provided card keys
        if "EDGE" in cards.keys():
            warn(f"Provided edge in cars will be overwritten by kwarg {edge}")
            cards.pop("EDGE")

        self._cards = cards
        self._radius = radius
        self._spectrum = spectrum
        self._name = name
        self._edge = edge
        self._feff_dict_set_kwargs = feff_dict_set_kwargs

        if self._edge == "L":
            warn(
                "Specified edge is 'L' and will be changed to L3 to be "
                "FEFF-compatible"
            )
            self._edge = "L3"
        elif "M" in self._edge:
            warn(
                "According to the FEFF9 documentaiton, M-edges are not well "
                "tested; proceed with caution!"
            )
        elif self._edge not in ["K", "L1", "L2", "L3"]:
            warn(
                f"Provided edge {self._edge} is not one of the standard "
                "choices K, L1, L2, L3 and it is unknown if FEFF supports it"
            )

    def get_FEFFDictSets(self, structure, absorbing_sites):
        """Constructs and returns a list of the
        :class:`pymatgen.io.feff.sets.FEFFDictSet` objects.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            The Pymatgen structure. Note that the ``absorbing_sites`` must
            correspond to the provided structure.
        absorbing_sites : list
            A list of int corresponding to absorbing sites.

        Returns
        -------
        list
            A list of :class:`pymatgen.io.feff.sets.FEFFDictSet` objects. One
            for each absorbing site.

        Raises
        ------
        ValueError
            If an invalid ``spectrum`` argument is provided.
        """

        if self._spectrum == "XANES":
            default_cards = copy(MPXANESSet.CONFIG)
        elif self._spectrum == "EXAFS":
            default_cards = copy(MPEXAFSSet.CONFIG)
        else:
            raise ValueError(f"Unknown spectrum type {self._spectrum}")

        default_cards.pop("EDGE")
        default_keys = default_cards.keys()
        keys_to_remove = list(set(default_keys) - set(self._cards.keys()))
        user_tag_settings = {"_del": keys_to_remove}
        user_tag_settings.update(self._cards)

        return [
            FEFFDictSet(
                site,
                structure,
                spectrum=self._spectrum,
                config_dict=default_cards,
                edge=self._edge,
                radius=self._radius,
                user_tag_settings=user_tag_settings,
                **self._feff_dict_set_kwargs,
            )
            for site in absorbing_sites
        ]

    def write(self, target_directory, **kwargs):
        """Writes the input files for the provided structure and sites. In the
        case of FEFF, if sites is None (usually indicating a global calculation
        such as a neutral potential electronic relaxation method in VASP), then
        write does nothing.

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the FEFF input files.
        **kwargs
            Must contain the ``structure`` key (the
            :class:`pymatgen.core.structure.Structure` of interest) and the
            ``sites`` key (a list of int, where each int corresponds to the
            site index of the site to write).

        Returns
        -------
        dict
            A dictionary containing the status and errors key. In the case of
            FEFF, there are no possible errors at this stage other than
            critical ones that would cause program termination, so the returned
            object is always
            ``{"pass": True, "errors": dict(), "paths": [...]}``.
        """

        structure = kwargs["structure_uc"]
        sites = kwargs["sites"]

        target_directory = Path(target_directory)
        target_directory.mkdir(exist_ok=True, parents=True)

        # Get the directory names
        species = [structure[site].specie.symbol for site in sites]
        names = [f"{site:03}_{specie}" for site, specie in zip(sites, species)]

        dict_sets = self.get_FEFFDictSets(structure, sites)

        paths = []
        for dict_set, name in zip(dict_sets, names):
            path = target_directory / Path(name)
            dict_set.write_input(path)
            paths.append(path)

        return {"pass": True, "errors": dict(), "paths": paths}

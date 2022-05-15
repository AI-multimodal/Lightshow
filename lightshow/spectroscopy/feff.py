from copy import copy

from monty.json import MSONable
from pymatgen.io.feff.sets import MPXANESSet, MPEXAFSSet, FEFFDictSet


class FEFFParameters(MSONable):
    """A one-stop-shop for all the different ways to modify input parameters
    for a FEFF calculation. This class is a lightweight wrapper for the
    ``FEFFDictSet``, containing some extra information and methods used for
    writing the appropriate input files. Like most Pymatgen objects, the
    class is serializable via e.g. ``feff_parameters.as_dict()``.

    Parameters
    ----------
    absorbing_sites : list
        A list of integers corresponding to absorbing sites. Note that these
        should be determined ahead of time.
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
    nkpts : int
        The number of k-points used in the Brillouin zone. Only used if
        FEFF is run in reciprocal space mode.
    """

    @property
    def absorbing_sites(self):
        return self._absorbing_sites

    @property
    def cards(self):
        return self._cards

    @property
    def edge(self):
        return self._edge

    @property
    def radius(self):
        return self._radius

    @property
    def nkpts(self):
        return self._nkpts

    def get_FEFFDictSets(self, structure, spectrum="XANES"):
        """Constructs and returns a list of the
        :class:`pymatgen.io.feff.sets.FEFFDictSet objects.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            The Pymatgen structure. Note that the ``absorbing_sites`` must
            correspond to the provided structure.
        spectrum : str, optional
            The type of spectroscopy to be run. This is used to set the default
            cards through Pymatgen. These defaults can be overridden by setting
            the ``cards`` argument when instantiating the class.

        Returns
        -------
        list
            A list of FEFFDictSet objects. One for each absorbing site.

        Raises
        ------
        ValueError
            If an invalid ``spectrum`` argument is provided.
        """

        if spectrum == "XANES":
            default_cards = copy(MPXANESSet.CONFIG)
        elif spectrum == "EXAFS":
            default_cards = copy(MPEXAFSSet.CONFIG)
        else:
            raise ValueError(f"Unknown spectrum type {spectrum}")

        default_cards.pop("EDGE")
        default_keys = default_cards.keys()
        keys_to_remove = list(set(default_keys) - set(self._cards.keys()))
        user_tag_settings = {"_del": keys_to_remove}
        user_tag_settings.update(self._cards)

        return [
            FEFFDictSet(
                site,
                structure,
                spectrum=spectrum,
                config_dict=default_cards,
                edge=self._edge,
                radius=self._radius,
                nkpts=self._nkpts,
                user_tag_settings=user_tag_settings,
            )
            for site in self._absorbing_sites
        ]

    def __init__(
        self, absorbing_sites, cards, edge="K", radius=9.0, nkpts=1000
    ):
        self._absorbing_sites = absorbing_sites
        self._cards = cards
        self._edge = edge
        self._radius = radius
        self._nkpts = nkpts

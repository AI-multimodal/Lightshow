from copy import copy
from pathlib import Path

from pymatgen.io.feff.sets import MPXANESSet, MPEXAFSSet, FEFFDictSet

from lightshow.parameters._base import _BaseParameters


class FEFFParameters(_BaseParameters):
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
    nkpts : int
        The number of k-points used in the Brillouin zone. Only used if
        FEFF is run in reciprocal space mode.
    spectrum : str
        The type of spectroscopy to be run. This is used to set the default
        cards through Pymatgen. These defaults can be overridden by setting
        the ``cards`` argument when instantiating the class. The type of
        calculation to run. Should likely be either ``"XANES"`` or ``"EXAFS"``.
    """

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

    @property
    def spectrum(self):
        return self._spectrum

    @property
    def calculation_name(self):
        """This is the name of the directory that will correspond to the type
        of calculation being run. In this case, FEFF-{spectrum}.

        Returns
        -------
        str
        """

        return f"FEFF-{self._spectrum}"

    def __init__(
        self,
        cards,
        edge="K",
        radius=9.0,
        nkpts=1000,
        spectrum="XANES",
    ):
        self._cards = cards
        self._edge = edge
        self._radius = radius
        self._nkpts = nkpts
        self._spectrum = spectrum

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
                nkpts=self._nkpts,
                user_tag_settings=user_tag_settings,
            )
            for site in absorbing_sites
        ]

    def validate(self, structure, sites):
        """Validates that the structure is compatible with a FEFF calculation.
        For FEFF, this method always returns True, as there is nothing specific
        to check.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            The Pymatgen structure of interest.

        Returns
        -------
        bool
            Always returns True.
        """

        return True

    def write(self, root, structure, sites):
        """Writes the input files for the provided structure and sites. In the
        case of FEFF, if sites is None (usually indicating a global calculation
        such as a neutral potential electronic relaxation method in VASP), then
        write does nothing.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            The Pymatgen structure of interest.
        sites : list
            A list of int, where each int corresponds to the site index of the
            site to write.
        """

        root = Path(root)
        root.mkdir(exist_ok=True, parents=True)

        # Get the directory names
        species = [structure[site].specie.symbol for site in sites]
        names = [f"{site:03}_{specie}" for site, specie in zip(sites, species)]

        dict_sets = self.get_FEFFDictSets(structure, sites)

        for dict_set, name in zip(dict_sets, names):
            path = root / Path(name)
            dict_set.write_input(path)

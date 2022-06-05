from pathlib import Path
import numpy as np

from monty.json import MSONable
from pymatgen.io.exciting import ExcitingInput

from lightshow.parameters._base import _BaseParameters


class EXCITINGParameters(MSONable, _BaseParameters):
    """A one-stop-shop for all the different ways to modify input parameters
    for an Exciting calculation. This class is a lightweight wrapper for the
    ``ExcitingInput``, containing some extra information and methods used for
    writing the appropriate input files. Like most Pymatgen objects, the
    class is serializable via e.g. ``feff_parameters.as_dict()``. !! TODO

    Parameters
    ----------
    cards : dict
        A dictionary of of the cards to be control the paramters in the
        EXCITING calculations.
        For example, one might wish to use something like

        .. code-block:: python
        cards = {
            'groundstate': {
                "xctype": "GGA_PBE",
                "nempty": "30",
                "rgkmax": "9.0",
                "do": "skip"
                },
            'xs': {
                "xstype": "BSE",
                "vkloff": "0.05 0.03 0.13",
                "nempty": "30",
                "gqmax": "4.0",
                "broad": "0.0327069",
                "tevout": "true",
                "tappinfo": "true",
                "energywindow": {"intv": "178.2 180.5",
                                 "points": "1000"},
                "screening": {"screentype": "full",
                              "nempty": "100"},
                "BSE": {"xas": "true",
                        "xasedge": "K",
                        "bsetype": "singlet",
                        "nstlxas": "1 20"}
                },

    """

    @property
    def name(self):
        return self._name

    def __init__(
        self,
        cards={
            "groundstate": {
                "xctype": "GGA_PBE",
                "nempty": "30",
                "rgkmax": "9.0",
                "do": "skip",
            },
            "xs": {
                "xstype": "BSE",
                "vkloff": "0.05 0.03 0.13",
                "nempty": "30",
                "gqmax": "4.0",
                "broad": "0.0327069",
                "tevout": "true",
                "tappinfo": "true",
                "energywindow": {"intv": "178.2 180.5", "points": "1000"},
                "screening": {"screentype": "full", "nempty": "100"},
                "BSE": {
                    "xas": "true",
                    "xasedge": "K",
                    "bsetype": "singlet",
                    "nstlxas": "1 20",
                },
            },
        },
        kpoints_method="custom",
        kpoints_method_kwargs={"cutoff": 32.0, "max_radii": 50.0},
        name="EXCITING",
    ):
        self._cards = cards
        # Method for determining the kmesh
        self._kpoints_method = kpoints_method
        self._kpoints_method_kwargs = kpoints_method_kwargs
        self._name = name

    @staticmethod
    def _getCondBands(volume, eRange):
        """Return a guess at the number of conduction bands that a given
        unit-cell volume needs to cover a given energy range (in Ryd).

        Parameters
        ----------
        volume : int
            The volume of the structure
        eRange : float
            The desired energy range

        Returns
        -------
        int :
            A number for the esimation number of bands
        """
        return round(0.256 * volume * (eRange ** (3 / 2)))

    @staticmethod
    def _getKmesh(structure, cutoff=32.0, max_radii=50.0):
        """Customizes the kmesh. For a kmesh sampling, e.g. [m, n, p],
        of a crystal "cell" is equivalent to generating a supercell with
        [m, n, p] the crystal cell. The corresponding radius is the largest
        radius of the sphere that can fit into this supercell. The radius can
        also be regarded as the inverse kmesh density. Along each reciprocal
        lattice, the kmesh densities are not guaranteed to be the same. The
        smallest kmesh density is chosen. This function uses the effective
        radius as a controlling factor to
        determining the kemsh.

        Parameters
        ----------
        supercell : pymatgen.core.structure.Structure
            Previously generated supercell (or standard unit cell).
        cutoff : float, optional
            Cutoff radius for constructing the kmesh. It will loop a look-up
            table for the effective radius (controlled by max_radii). The kmesh
            with radius right above the cutoff will be chosen. Default is 32
            Angstroms (60 Bohr).
        max_radii : float, optional
            Maximum radius used for constructing the lookup table.

        Returns
        -------
        tuple
        """

        klist = dict()
        rlatt = np.array(structure.lattice.reciprocal_lattice.abc)
        # TODO: why 10, 0.2?
        for xx in np.arange(0, 10, 0.2):
            div = np.floor(xx * rlatt) + 1
            divlatt = 2.0 * np.pi / rlatt * div

            radi = min(divlatt)

            if radi > max_radii:
                break
            else:
                div = tuple(div.astype(int))
                if div not in klist:
                    klist[div] = radi

        for key, value in klist.items():
            if value > cutoff:
                k = key
                break

        return k

    def write(self, target_directory, **kwargs):
        """Writes the input files for the provided structure and sites. In the
        case of FEFF, if sites is None (usually indicating a global calculation
        such as a neutral potential electronic relaxation method in VASP), then
        write does nothing. # TODO

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the FEFF input files.
        **kwargs
            Must contain the ``structure_uc`` key (the
            :class:`pymatgen.core.structure.Structure` of interest) and the
            ``sites`` key (a list of int, where each int corresponds to the
            site index of the site to write).

        Returns
        -------
        dict
            A dictionary containing the status and errors key. In the case of
            EXCITING, there are no possible errors at this stage other than
            critical ones that would cause program termination, so the returned
            object is always ``{"pass": True, "errors": dict()}``.
        """

        structure = kwargs["structure_uc"]
        sites = kwargs["sites"]

        target_directory = Path(target_directory)
        target_directory.mkdir(exist_ok=True, parents=True)

        # Get the directory names
        species = [structure[site].specie.symbol for site in sites]

        excitinginput = ExcitingInput(structure)
        # Estimate number of band
        nbands = self._getCondBands(structure.lattice.volume, 2.25)
        self._cards["xs"]["BSE"]["nstlxas"] = f"1 {nbands}"
        # Estimate number of kpoints
        kmesh = self._getKmesh(structure, cutoff=16.0, max_radii=50.0)
        self._cards["groundstate"][
            "ngridk"
        ] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        self._cards["xs"]["ngridk"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        self._cards["xs"]["ngridq"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"

        # Determine XAS species
        species = [structure[site].specie.symbol for site in sites]

        i = 0
        for specie in sorted(structure.types_of_species, key=lambda el: el.X):
            i = i + 1
            if specie.symbol == species[0]:
                self._cards["xs"]["BSE"]["xasspecies"] = str(i)

        for site, specie in zip(sites, species):
            path = target_directory / Path(f"{site:03}_{specie}")
            path.mkdir(exist_ok=True, parents=True)

            filepath_xas = path / "input_xas.xml"
            self._cards["xs"]["BSE"]["xasatom"] = str(site + 1)
            excitinginput.write_file(
                "primitive", filepath_xas, bandstr=False, **self._cards
                )

        return {"pass": True, "errors": dict()}

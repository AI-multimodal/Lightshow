from pathlib import Path

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
                    "energywindow": {
                        "intv": "178.2 180.5",
                        "points": "1000"
                    },
                    "screening": {
                        "screentype": "full",
                        "nempty": "100"
                    },
                    "BSE": {
                        "xas": "true",
                        "xasedge": "K",
                        "bsetype": "singlet",
                        "nstlxas": "1 20"
                    }
                }
            }

    """

    @property
    def name(self):
        return self._name

    def __init__(
        self,
        cards={
            "structure": {"speciespath": "./", "autormt": "true"},
            "groundstate": {
                "xctype": "GGA_PBE",
                "nempty": "200",
                "rgkmax": "9.0",
                "do": "fromscratch",
            },
            "xs": {
                "xstype": "BSE",
                "vkloff": "0.05 0.03 0.13",
                "nempty": "150",
                "gqmax": "4.0",
                "broad": "0.0327069",
                "tevout": "true",
                "tappinfo": "true",
                "energywindow": {"intv": "178.2 180.5", "points": "1000"},
                "screening": {"screentype": "full", "nempty": "150"},
                "BSE": {
                    "xas": "true",
                    "xasedge": "K",
                    "bsetype": "singlet",
                    "nstlxas": "1 20",
                    "distribute": "true",
                    "eecs": "1000",
                },
                "qpointset": {"qpoint": {"text()": "0.0 0.0 0.0"}},
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

            filepath_xas = path / "input.xml"
            self._cards["xs"]["BSE"]["xasatom"] = str(site + 1)
            excitinginput.write_file(
                "primitive", filepath_xas, bandstr=False, **self._cards
            )

        return {"pass": True, "errors": dict()}
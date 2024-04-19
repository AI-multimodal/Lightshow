import xml.etree.ElementTree as ET
from pathlib import Path
from warnings import warn

from monty.json import MSONable
from pymatgen.io.exciting import ExcitingInput

from lightshow.common.kpoints import GenericEstimatorKpoints
from lightshow.common.nbands import UnitCellVolumeEstimate
from lightshow.parameters._base import _BaseParameters
from lightshow.utils.environ_utils import get_SPECIES_DIRECTORY_from_environ

EXCITING_DEFAULT_GQMAX = 4.0

EXCITING_DEFAULT_CARDS = {
    "structure": {"speciespath": "./", "autormt": "true"},
    "groundstate": {
        "xctype": "GGA_PBE",
        "nempty": "200",
        "rgkmax": "9.0",
        "do": "fromscratch",
        "gmaxvr": "25",
        "lmaxmat": "10",
        "lmaxvr": "10",
        "lmaxapw": "10",
    },
    "xs": {
        "xstype": "BSE",
        "vkloff": "0.05 0.03 0.13",
        "nempty": "150",
        "gqmax": str(EXCITING_DEFAULT_GQMAX),
        "broad": "0.0327069",
        "tevout": "true",
        "tappinfo": "true",
        "energywindow": {"intv": "178.2 180.5", "points": "1000"},
        "screening": {"screentype": "full", "nempty": "150"},
        "BSE": {
            "xas": "true",
            "bsetype": "singlet",
            "nstlxas": "1 20",
            "distribute": "true",
            "eecs": "1000",
        },
        "qpointset": {"qpoint": {"text()": "0.0 0.0 0.0"}},
    },
}


class EXCITINGParameters(MSONable, _BaseParameters):
    """A one-stop-shop for all the different ways to modify input parameters
    for an Exciting calculation. This class is a lightweight wrapper for the
    ``ExcitingInput``, containing some extra information and methods used for
    writing the appropriate input files. Like most Pymatgen objects, the
    class is serializable via e.g. ``exciting_parameters.as_dict()``.

    Parameters
    ----------
    cards : dict
        A dictionary of of the cards to be control the paramters in the
        EXCITING calculations. It contains all the default setting in
        this workflow.

        .. code-block:: python

            cards = {
                "structure": {"speciespath": "./", "autormt": "true"},
                "groundstate": {
                    "xctype": "GGA_PBE",
                    "nempty": "200",
                    "rgkmax": "9.0",
                    "do": "fromscratch",
                    "gmaxvr": "25",
                    "lmaxmat": "10",
                    "lmaxvr": "10",
                    "lmaxapw": "10",
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
                }
            }
    edge : str
        The XAS edge for the calculation
    species_directory : str
        A string contains the absolute path for species files. Usually, this
        folder comes with the exciting package, which can be obtained at
        http://exciting.wikidot.com/
    kpoints : lightshow.common.kpoints._BaseKpointsMethod
        The method for constructing the kpoints file from the structure. Should
        be a class with a ``__call__`` method defined. This method should take
        the structure as input and return a tuple corresponding to the kpoints
        density along each axis.
    nbands : lightshow.common.nbands._BaseNbandsMethod
        The method for determining the number of valence bands from the
        structure. Should be a class with a ``__call__`` method defined. This
        method should take the structure as input and return an integer: the
        number of valence bands to use in the calculation.
    gqmax : float
        :math:`|G+q|` cutoff of the plane wave expansion
    name : str
        The name of the calculation. Should likely always be ``"EXCITING"``.
    """

    @property
    def name(self):
        return self._name

    def __init__(
        self,
        cards=EXCITING_DEFAULT_CARDS,
        species_directory=None,
        plan=[
            "xsgeneigvec",
            "writepmatxs",
            "scrgeneigvec",
            "scrwritepmat",
            "screen",
            "scrcoulint",
            "exccoulint",
            "bse",
        ],
        edge="K",
        kpoints=GenericEstimatorKpoints(cutoff=16.0, max_radii=50.0),
        nbands=UnitCellVolumeEstimate(e_range=30.0),
        gqmax=EXCITING_DEFAULT_GQMAX,
        name="EXCITING",
    ):
        # Default cards
        self._cards = cards

        # species directory
        if species_directory is None:
            species_directory = get_SPECIES_DIRECTORY_from_environ()
        if species_directory is None:
            warn(
                "species_directory not set, and SPECIES_DIRECTORY not in "
                "the current environment variables. The current/working "
                "folder will be used as default. Please make sure you copy "
                "all the corresponding species file into the working folder, "
                "e.g. where the input.xml file is generated."
            )
            species_directory = "./"
        self._species_directory = species_directory
        self._cards["structure"]["speciespath"] = self._species_directory
        # Update plan
        self._plan = plan

        # Set XAS edge
        self._edge = edge
        self._cards["xs"]["BSE"]["xasedge"] = self._edge
        # Method for determining the kmesh
        self._kpoints = kpoints
        # Method for determining number of bands
        self._nbands = nbands

        # Set gqmax
        self._gqmax = str(gqmax)
        self._cards["xs"]["gqmax"] = self._gqmax

        # Method for determining the kmesh
        self._name = name

    def write(self, target_directory, **kwargs):
        """Writes the input files for the provided structures and sites.

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

        excitinginput = ExcitingInput(structure)
        # Estimate number of band
        nbands = self._nbands(structure)
        self._cards["xs"]["BSE"]["nstlxas"] = f"1 {nbands}"
        # Estimate number of kpoints
        kmesh = self._kpoints(structure)
        self._cards["groundstate"]["ngridk"] = (
            f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        )
        self._cards["xs"]["ngridk"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        self._cards["xs"]["ngridq"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        # Determine XAS species
        species = [structure[site].specie.symbol for site in sites]

        for i, specie in enumerate(
            sorted(structure.types_of_species, key=lambda el: el.X)
        ):
            if specie.symbol == species[0]:
                self._cards["xs"]["BSE"]["xasspecies"] = str(i + 1)

        for site, specie in zip(sites, species):
            path = target_directory / f"{site:03}_{specie}"
            path.mkdir(exist_ok=True, parents=True)

            filepath_xas = path / "input.xml"
            self._cards["xs"]["BSE"]["xasatom"] = str(site + 1)
            # add the plans using the tree
            root = excitinginput.write_etree(
                "unchanged", bandstr=False, **self._cards
            )
            tree = ET.ElementTree(root)
            xs_loc = root.find("xs")
            plan = ET.Element("plan")
            xs_loc.insert(-1, plan)
            ploc = xs_loc.find("plan")
            for task in self._plan:
                ET.SubElement(ploc, "doonly", task=task)

            excitinginput._indent(root)
            tree.write(filepath_xas)

        return {"pass": True, "errors": dict()}

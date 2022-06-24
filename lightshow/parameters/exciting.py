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
    gs_cards : dict
        A dictionary of of the cards to be control the paramters in the
        EXCITING calculations: ground state part. Default is empty, which means
        it will use all the default setting in this workflow.
        For example, one might wish to use something like

        .. code-block:: python

            gs_cards = {
                "xctype": "GGA_PBE",
                "nempty": "30",
                "rgkmax": "9.0",
                "do": "skip"
                }
    xs_cards : dict
        A dictionary of of the cards to be control the paramters in the
        EXCITING calculations: xs part. Default is empty, which means
        it will use all the default setting in this workflow.
        For example, one might wish to use something like

        .. code-block:: python

            xs_cards = {
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
    species_path : str
        A string contains the absolute path for species files.
    kpoints_method : str
        Methods for determining the kmesh. Currently, only "custom" is supported.
    kpoints_method_kwargs : dict
        Arguments to pass to the classmethod used to construct the kpoints.
    nbands_estimator : str
        Methods for determining the kmesh. If 'heg', homogeneous electron gas
        model is used to estimate the number of condcution band. If the input is
        a number, it will be used as the number of bands (posive number stands for
        the number of all bands; negative number stands for the number of conduction
        bands).
    nbands_estimator_kwargs: dict
        Used only when the nbands_estimator == 'heg'. Optional arguments to pass to
        the classmethod to estimate the number of conduction bands.
    """

    @property
    def name(self):
        return self._name

    def __init__(
        self,
        gs_cards={},
        xs_cards={},
        species_path="./",
        nbands_estimator="heg",
        nbands_estimator_kwargs={"eRange": 2.25},  # unit in Ryd
        kpoints_method="custom",
        kpoints_method_kwargs={"cutoff": 32.0, "max_radii": 50.0},
        name="EXCITING",
    ):
        # Default cards
        self._cards = {
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
            },
        }
        # Update speciespath
        self._species_path = species_path
        self._cards["structure"]["speciespath"] = self._species_path
        # User input parameters
        self._gs_cards = gs_cards
        self._xs_cards = xs_cards
        # Method for determining the nbands
        self._nbands_estimator = nbands_estimator
        self._nbands_estimator_kwargs = nbands_estimator_kwargs  # unit in Ryd
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
        if self._nbands_estimator == "heg":
            eRange = self._nbands_estimator_kwargs["eRange"]
            nbands = self._getCondBands(structure.lattice.volume, eRange)
            self._cards["xs"]["BSE"]["nstlxas"] = f"1 {nbands}"
        else:
            try:
                nbands = int(self._nbands_estimator)
                if nbands <= 0:
                    raise ValueError(
                        "the input of nbands_estimator is not supported"
                    )

                self._cards["xs"]["BSE"]["nstlxas"] = f"1 {nbands}"
            except Exception:
                raise ValueError(
                    "the input of nbands_estimator is not supported"
                )
        # Estimate number of kpoints
        if self._kpoints_method == "custom":
            cutoff = self._kpoints_method_kwargs["cutoff"]
            max_radii = self._kpoints_method_kwargs["max_radii"]
            kmesh = self._getKmesh(
                structure, cutoff=cutoff, max_radii=max_radii
            )
            self._cards["groundstate"][
                "ngridk"
            ] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
            self._cards["xs"]["ngridk"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
            self._cards["xs"]["ngridq"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        else:
            raise ValueError("method for obtaining kmesh not supported")
        # Update user's setting
        for (
            key,
            val,
        ) in self._gs_cards.items():  # only have two layers, which if more?
            if isinstance(val, dict):
                for sub_key, sub_val in val.items():
                    self._cards["groundstate"][key][sub_key] = sub_val
            else:
                self._cards["groundstate"][key] = val
        for key, val in self._xs_cards.items():
            if isinstance(val, dict):
                for sub_key, sub_val in val.items():
                    self._cards["xs"][key][sub_key] = sub_val
            else:
                self._cards["xs"][key] = val
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

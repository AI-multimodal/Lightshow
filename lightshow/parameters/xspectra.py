from pathlib import Path
import os
import json
import bz2
import base64
from collections import OrderedDict

from monty.json import MSONable
from pymatgen.io.pwscf import PWInput

from lightshow.parameters._base import _BaseParameters
from lightshow.common.kpoints import GenericEstimatorKpoints
import lightshow


XSPECTRA_DEFAULT_CARDS = {
    "QE": {
        "control": {"restart_mode": "from_scratch", "wf_collect": ".true."},
        "electrons": {"conv_thr": 1e-08, "mixing_beta": 0.4},
        "system": {
            "degauss": 0.002,
            "ecutrho": 320,
            "ecutwfc": 40,
            "nspin": 1,
            "occupations": "smearing",
            "smearing": "gauss",
        },
    },
    "XS": {
        "cut_occ": {"cut_desmooth": 0.3},
        "input_xspectra": {
            "edge": "K",
            "outdir": "../",
            "prefix": "pwscf",
            "xcheck_conv": 200,
            "xerror": 0.01,
            "xniter": 5000,
            "xcoordcrys": ".false.",
        },
        "kpts": {"kpts": "2 2 2", "shift": "0 0 0"},
        "plot": {
            "cut_occ_states": ".true.",
            "terminator": ".true.",
            "xemax": 70,
            "xemin": -15.0,
            "xnepoint": 400,
        },
        "plotFalse": {
            "gamma_energy(1)": 7,
            "gamma_energy(2)": 23,
            "gamma_mode": "variable",
            "gamma_value(1)": 0.89,
            "gamma_value(2)": 2.1,
        },
        "plotTrue": {"gamma_mode": "constant", "xgamma": 0.05},
    },
    "XS_controls": {
        "Rmin": "9.0",
        "psp_json": "SSSP_precision",
    },
}


class XSpectraParameters(MSONable, _BaseParameters):
    """A one-stop-shop for all the different ways to modify input parameters
    for an XSpectra calculation. !! TODO

    Parameters
    ----------
    cards : dict
        A dictionary of of the cards to be control the parameters in the
        OCEAN calculations. For example, one might wish to use something like

        .. code-block:: python

            cards = {
                "QE": {
                    "control": {
                        "restart_mode": "from_scratch",
                        "wf_collect": ".true."
                    },
                    "electrons": {
                        "conv_thr": 1e-08,
                        "mixing_beta": 0.4
                    },
                    "ions": {},
                    "system": {
                        "degauss": 0.002,
                        "ecutrho": 320,
                        "ecutwfc": 40,
                        "nspin": 1,
                        "occupations": "smearing",
                        "smearing": "gauss"
                    }
                },
                "XS": {
                    "cut_occ": {
                        "cut_desmooth": 0.3
                    },
                    "input_xspectra": {
                        "edge": "K",
                        "outdir": "../",
                        "prefix": "pwscf",
                        "!wf_collect": ".true.",
                        "xcheck_conv": 200,
                        "xerror": 0.01,
                        "xniter": 5000,
                        "xcoordcrys": ".false."
                    },
                    "kpts": {
                        "kpts": "2 2 2",
                        "shift": "0 0 0"
                    },
                    "plot": {
                        "cut_occ_states": ".true.",
                        "terminator": ".true.",
                        "xemax": 70,
                        "xemin": -15.0,
                        "xnepoint": 400
                    },
                    "plotFalse": {
                        "gamma_energy(1)": 7,
                        "gamma_energy(2)": 23,
                        "gamma_mode": "variable",
                        "gamma_value(1)": 0.89,
                        "gamma_value(2)": 2.1
                    },
                    "plotTrue": {
                        "gamma_mode": "constant",
                        "xgamma": 0.05
                    }
                },
                "XS_controls": {
                    "element": "Ti",
                    "edge": "K",
                    "kden": "-1",
                    "psp": {
                        "Ti+": "Ti.fch.upf",
                    },
                    "Rmin": "9.0",
                    "scf_kden": "-1",
                    "core_psp_json": "FCH1",
                    "psp_json": "SSSP_precision"
                }
            }
    kpoints : lightshow.common.kpoints._BaseKpointsMethod
        The method for constructing he kpoints file from the structure. Should
        be a class with a ``__call__`` method defined. This method should take
        the structure as input and return a tuple corresponding to the kpoints
        density along each axis.
    nbands : lightshow.common.nbands._BaseNbandsMethod
        The method for determining the number of valence bands from the
        structure. Should be a class with a ``__call__`` method defined. This
        method should take the structure as input and return an integer: the
        number of valence bands to use in the calculation.

    """

    @property
    def name(self):
        return self._name

    @property
    def cards(self):
        return self._cards

    def __init__(
        self,
        cards=XSPECTRA_DEFAULT_CARDS,
        kpoints=GenericEstimatorKpoints(cutoff=16.0, max_radii=50.0),
        defaultConvPerAtom=1e-10,
        edge="K",
        name="XSpectra",
    ):
        self._cards = cards
        # Method for determining the kmesh
        self._kpoints = kpoints
        self._defaultConvPerAtom = defaultConvPerAtom
        self._edge = edge
        self._name = name

    @staticmethod
    def _unpackPsps(
        ecutwfc,
        ecutrho,
        pspDatabaseRoot,
        DatabaseDir,
        symbols,
        folder,
        needWfn=False,
    ):
        psp = {}
        pspDatabaseName = pspDatabaseRoot + ".json"
        sssp_fn = os.path.join(DatabaseDir, pspDatabaseName)
        with open(sssp_fn, "r") as pspDatabaseFile:
            pspDatabase = json.load(pspDatabaseFile)
        minSymbols = set(symbols)
        for symbol in minSymbols:
            # print(symbol)
            # print(pspDatabase[symbol]["filename"])
            psp[symbol] = pspDatabase[symbol]["filename"]
            if ecutwfc < pspDatabase[symbol]["cutoff"]:
                ecutwfc = pspDatabase[symbol]["cutoff"]
            if ecutrho < pspDatabase[symbol]["rho_cutoff"]:
                ecutrho = pspDatabase[symbol]["rho_cutoff"]
        #        if xsJSON['QE']['system']['ecutwfc'] < pspDatabase[ symbol ]['cutoff']:
        #            xsJSON['QE']['system']['ecutwfc'] = pspDatabase[ symbol ]['cutoff']
        #        if xsJSON['QE']['system']['ecutrho'] < pspDatabase[ symbol ]['rho_cutoff']:
        #            xsJSON['QE']['system']['ecutrho'] = pspDatabase[ symbol ]['rho_cutoff']

        pspDatabaseName = pspDatabaseRoot + "_pseudos.json"
        sssp_fn = Path(DatabaseDir) / Path(pspDatabaseName)
        with open(sssp_fn, "r") as p:
            pspJSON = json.load(p)
        for symbol in minSymbols:
            fileName = psp[symbol]
            pspString = bz2.decompress(base64.b64decode(pspJSON[fileName]))
            # print("Expected hash:  " + pspDatabase[symbol]["md5"])
            # print("Resultant hash: " + hashlib.md5(pspString).hexdigest())
            with open(folder / fileName, "w") as f:
                f.write(pspString.decode("utf-8"))

        if needWfn:
            for symbol in minSymbols:
                if "wfc" not in pspDatabase[symbol]:
                    print(
                        "WFC not stored corectly in "
                        + pspDatabaseRoot
                        + " for element "
                        + symbol
                    )
                    return False
                fileName = pspDatabase[symbol]["wfc"]
                pspString = bz2.decompress(base64.b64decode(pspJSON[fileName]))
                # print("Expected hash:  " + pspDatabase[symbol]["wfc_md5"])
                # print("Resultant hash: " + hashlib.md5(pspString).hexdigest())
                element = symbol.split("+")[0]
                with open(folder / f"Core_{element}.wfc", "w") as f:
                    f.write(pspString.decode("utf-8"))

        return psp, ecutwfc, ecutrho

    @staticmethod
    def _write_xspectra_in(
        mode, iabs, dirs, xkvec, element, XSparams: dict, plot=False
    ):
        """construct input file for XSpectra calculation

        Parameters
        ----------
        mode : str, mandatory
            "dipole" or "quadrupole"
        iabs : int, mandatory
            the index of the absorbing element in scf calculation
        dirs : str, mandatory
            TODO
        xkvec : tuple, mandatory
            TODO
        XSparams : dict, mandatory
            paramers parsed to XSpectra calculation
        plot : boolen, optional
            controls xonly_plot

        Returns
        -------
        string of the XSpectra input file
        """
        inp = [
            "&input_xspectra",
            "    calculation = 'xanes_%s'" % mode,
            "    edge = '" + XSparams["input_xspectra"]["edge"] + "'",
            "    prefix = 'pwscf'",
            "    outdir = '../'",
            "    xniter = " + str(XSparams["input_xspectra"]["xniter"]),
            "    xiabs = %d" % iabs,
            "    xerror = " + str(XSparams["input_xspectra"]["xerror"]),
            "    xcoordcrys = '"
            + XSparams["input_xspectra"]["xcoordcrys"]
            + "'",
            "    xcheck_conv = "
            + str(XSparams["input_xspectra"]["xcheck_conv"]),
            "    xepsilon(1) = %d" % dirs[0],
            "    xepsilon(2) = %d" % dirs[1],
            "    xepsilon(3) = %d" % dirs[2],
        ]

        if mode == "quadrupole":
            inp += [
                "    xkvec(1) = %.10f" % xkvec[0],
                "    xkvec(2) = %.10f" % xkvec[1],
                "    xkvec(3) = %.10f" % xkvec[2],
            ]

        if plot:
            inp += ["    xonly_plot = .true."]

        inp += [
            "/",
            "&plot",
            "    xnepoint = " + str(XSparams["plot"]["xnepoint"]),
            "    xemin = " + str(XSparams["plot"]["xemin"]),
            "    xemax = " + str(XSparams["plot"]["xemax"]),
            "    terminator = " + XSparams["plot"]["terminator"],
            "    cut_occ_states = " + XSparams["plot"]["cut_occ_states"],
        ]
        if plot:
            inp += [
                "    xgamma = " + str(XSparams["plotTrue"]["xgamma"]),
                "    gamma_mode = '"
                + XSparams["plotTrue"]["gamma_mode"]
                + "'",
                "/",
            ]
        else:
            # use very small smearing value: 0.01 eV
            # Table IIA in Campbell and Papp (2001) https://doi.org/10.1006/adnd.2000.0848
            inp += ["    gamma_mode = 'constant'", "    xgamma = 0.01 ", "/"]

        inp += [
            "&pseudos",
            f"    filecore = '../../Core_{element}.wfc'",
            # "    r_paw(1) = 1.79",  # hard-coded to Ti w/ core-hole
            "/",
            "&cut_occ",
            "    cut_desmooth = " + str(XSparams["cut_occ"]["cut_desmooth"]),
            "/",
            XSparams["kpts"]["kpts"] + " " + XSparams["kpts"]["shift"],
        ]
        return "\n".join(inp) + "\n"

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

        structure = kwargs["structure_sc"]
        sites = kwargs["sites"]
        index_mapping = kwargs["index_mapping"]

        target_directory = Path(target_directory)
        target_directory.mkdir(exist_ok=True, parents=True)
        symbols = [spec.symbol for spec in structure.species]
        # Obtain absorbing atom
        species = [
            structure[index_mapping[site]].specie.symbol for site in sites
        ]
        element = species[0]
        self._cards["XS_controls"]["element"] = element
        self._cards["XS_controls"]["edge"] = self._edge
        symTarg = self._cards["XS_controls"]["element"]
        # Estimate number of kpoints
        if float(self._cards["XS_controls"]["Rmin"]) >= 9:
            # use Gamma point for ground state calculations (es.in and gs.in)
            kpoints_scf = [1, 1, 1]
        else:
            kpoints_scf = self._kpoints(structure)

        kpoints_xas = self._kpoints(structure)

        self._cards["XS"]["kpts"][
            "kpts"
        ] = f"{kpoints_xas[0]} {kpoints_xas[1]} {kpoints_xas[2]}"
        # Determine the SCF? convergence threshold
        self._cards["QE"]["electrons"][
            "conv_thr"
        ] = self._defaultConvPerAtom * len(structure)
        # Get the psp data ready for the GS calculations; similar to SCF (neutral) calculations in VASP
        module_path = Path(lightshow.parameters.__path__[0])
        pspDatabaseRoot = self._cards["XS_controls"]["psp_json"]
        DatabaseDir = module_path / "pseudos" / "data"
        ecutwfc = self._cards["QE"]["system"]["ecutwfc"]
        ecutrho = self._cards["QE"]["system"]["ecutrho"]
        psp, ecutwfc, ecutrho = self._unpackPsps(
            ecutwfc,
            ecutrho,
            pspDatabaseRoot,
            DatabaseDir,
            symbols,
            target_directory,
        )

        self._cards["QE"]["system"]["ecutwfc"] = ecutwfc
        self._cards["QE"]["system"]["ecutrho"] = ecutrho
        self._cards["QE"]["control"]["pseudo_dir"] = "../"

        path = target_directory / "GS"
        path.mkdir(exist_ok=True, parents=True)
        gs_in = PWInput(
            structure,
            pseudo=psp,
            control=self._cards["QE"]["control"],
            system=self._cards["QE"]["system"],
            electrons=self._cards["QE"]["electrons"],
            kpoints_grid=kpoints_scf,
        )
        gs_in.write_file(path / "gs.in")
        # Get the psp data read for ES calculations
        # pspDatabaseRoot = self._cards["XS_controls"]["core_psp_json"]
        # try:
        #    psp2, ecutwfc, ecutrho = self._unpackPsps(
        #        ecutwfc,
        #        ecutrho,
        #        pspDatabaseRoot,
        #        DatabaseDir,
        #        [self._cards["XS_controls"]["element"]],
        #        target_directory,
        #        needWfn=True,
        #    )
        # except KeyError:
        # throw a warning here
        # warn("take care of the core-hole psp for absorber by yourself")

        # give a name to the pseudo potential
        # set relatively large cutoff as default
        # psp2 = {element: f"{element}.fch.upf"}
        # ecutwfc = 100
        # ecutrho = 800
        # for i in psp2:
        psp[f"{element}+"] = f"{element}.fch.upf"  # psp2[i]
        # Determine iabs
        psp = OrderedDict(psp)
        for i, j in enumerate(psp.keys()):
            if j == symTarg + "+":
                iabs = i + 1
        for site, specie in zip(sites, species):
            path = target_directory / Path(f"{site:03}_{specie}")
            path.mkdir(exist_ok=True, parents=True)

            structure[index_mapping[site]] = element + "+"
            self._cards["QE"]["control"]["pseudo_dir"] = "../"
            es_in = PWInput(
                structure,
                pseudo=psp,
                control=self._cards["QE"]["control"],
                system=self._cards["QE"]["system"],
                electrons=self._cards["QE"]["electrons"],
                kpoints_grid=kpoints_scf,
            )
            es_in.write_file(path / "es.in")
            structure[index_mapping[site]] = element

            # Deal with the dipole case only
            # notice I put the photonSymm in the folder, which is created by John
            photons = list()
            photons.append({"dipole": [1, 0, 0, 1]})
            photons.append({"dipole": [0, 1, 0, 1]})
            photons.append({"dipole": [0, 0, 1, 1]})

            totalweight = 0
            for photon in photons:
                totalweight += photon["dipole"][3]

            photonCount = 0
            for photon in photons:
                photonCount += 1
                dir1 = photon["dipole"][0:3]
                dir2 = dir1
                weight = photon["dipole"][3] / totalweight
                mode = "dipole"

                xanesfolder = path / f"{mode}{photonCount}"
                xanesfolder.mkdir(parents=True, exist_ok=True)
                with open(xanesfolder / "xanes.in", "w") as f:
                    f.write(
                        self._write_xspectra_in(
                            mode,
                            iabs,
                            dir1,
                            dir2,
                            self._cards["XS_controls"]["element"],
                            self._cards["XS"],
                        )
                    )

                with open(xanesfolder / "weight.txt", "w") as f:
                    f.write(str(weight) + "\n")

        return {"pass": True, "errors": dict()}

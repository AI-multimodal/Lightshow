from copy import copy
from pathlib import Path
import numpy as np

from monty.json import MSONable
from pymatgen.core import Element
from ase.units import Bohr

from lightshow.parameters._base import _BaseParameters


OCEAN_DEFAULT_CARDS = {}


class OCEANParameters(MSONable, _BaseParameters):
    """A one-stop-shop for all the different ways to modify input parameters
    for an OCEAN calculation. !! TODO

    Parameters
    ----------
    input_cards : dict
        A dictionary of of the cards to be control the parameters in the
        OCEAN calculations. Default is empty, which means that users will
        use most default parameters stored at self._cards. If users want to
        change some default parameters, for example, one might wish to abinit
        rather than qe as DFT package, also they want to use the energy range
        for DFT to be 50 eV, they can something like

        .. code-block:: python

            cards = {
                "dft": "qe",
                "dft_energy_range": 50,
            }
    kpoints_method : str
        Methods for determining the kmesh. Currently, only "custom" is supported.
    kpoints_method_kwargs : dict
        Arguments to pass to the classmethod used to construct the kpoints.
    defaultConvPerAtom: float, Default: 1e-10 Ryd/atom
        Enegy convergence threshold.
    bandgap: float, Default: None
        User can provide the band gap of the material they are interested in
        to determine the `diemac`. If not provided and the strucutre is from
        pymatgen, the code will use the data (either band gap or diel) stored
        in the database to get `diemac`; if not provided and the structure
        is from a file, the default value of `diemac` in the cards will be used.
        If `diel` is also not None, `diemac` is determined using `diel`
    diel: float, Default: None
        User can provide the diel of the material they are interested in
        to determine the `diemac`. If not provided and the strucutre is from
        pymatgen, the code will use the the data (either band gap or diel) stored
        in the database to get `diemac`; if not provided and the structure is
        from a file, the default value of `diemac` in the cards will be used.
    edge: str, Default: 'K'
        Type of edges
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

    @property
    def cards(self):
        return self._cards

    def __init__(
        self,
        cards=OCEAN_DEFAULT_CARDS,
        bandgap=None,
        diel=None,
        nbands_estimator="heg",
        nbands_estimator_kwargs={"eRange": 2.25},  # unit in Ryd
        kpoints_method="custom",
        kpoints_method_kwargs={"cutoff": 16.0, "max_radii": 50.0},
        defaultConvPerAtom=1e-10,
        edge="K",
        name="OCEAN",
    ):

        # Default cards for ocean
        self._cards = cards

        # User modified cards
        self._bandgap = bandgap
        self._diel = diel
        # Method for determining number of bands
        self._nbands_estimator = nbands_estimator
        self._nbands_estimator_kwargs = nbands_estimator_kwargs
        # Method for determining the kmesh
        self._kpoints_method = kpoints_method
        self._kpoints_method_kwargs = kpoints_method_kwargs
        self._defaultConvPerAtom = defaultConvPerAtom

        # Handle the edge parsing here.
        if edge not in self._edge_map.keys():
            raise ValueError(
                f"Provided edge {edge} is not a valid choice for OCEAN "
                "calculations. Edges should be chosen from the list: "
                f"{list(self._edge_map.keys())}"
            )
        self._edge = self._edge_map[edge]
        self._name = name

    @property
    def _edge_map(self):
        """Mapping between letter (indicating the edge) and number (input to
        OCEAN) defining the XAS edge."""

        return {
            "K": "1 0",
            "L": "2 1",
            "L1": "2 0",
            "L2": "2 1",
            "L3": "2 1",
            "L23": "2 1",
            "M": "3 2",
            "M1": "3 0",
            "M2": "3 1",
            "M3": "3 1",
            "M23": "3 1",
            "M4": "3 2",
            "M5": "3 2",
            "M45": "3 2",
            "N": "4 3",
            "N1": "4 0",
            "N2": "4 1",
            "N3": "4 1",
            "N23": "4 1",
            "N4": "4 2",
            "N5": "4 2",
            "N45": "4 2",
            "N6": "4 3",
            "N7": "4 3",
            "N67": "4 3",
        }

    @staticmethod
    def _oceanKptSampling(cell, kpt):
        # John's method for getting the explicit screen.nkpt
        # Should normalize each cell dim by Bohr, or we can just do it the once for volume
        v = abs(np.dot(np.cross(cell[0], cell[1]), cell[2])) / Bohr
        b1 = 2 * np.pi * np.linalg.norm(np.cross(cell[1], cell[2])) / v
        b2 = 2 * np.pi * np.linalg.norm(np.cross(cell[0], cell[2])) / v
        b3 = 2 * np.pi * np.linalg.norm(np.cross(cell[0], cell[1])) / v

        k1 = int(-kpt * b1) + 1
        k2 = int(-kpt * b2) + 1
        k3 = int(-kpt * b3) + 1

        return "{:d} {:d} {:d}".format(k1, k2, k3)

    @staticmethod
    def _write_ocean_in(path, structure, input_data: dict):

        fd = open(path, "w")

        input_data_str = []
        for key in input_data:
            input_data_str.append(
                str(key) + " { " + str(input_data[key]) + " }\n"
            )

        fd.write("".join(input_data_str))
        species = sorted(set(structure.atomic_numbers))

        fd.write("znucl {{ {} }}\n".format(" ".join(str(Z) for Z in species)))
        fd.write("typat")
        fd.write("{\n")
        types = []
        for Z in structure.atomic_numbers:
            for n, Zs in enumerate(species):
                if Z == Zs:
                    types.append(n + 1)
        n_entries_int = 20  # integer entries per line
        for n, type in enumerate(types):
            fd.write(" %d" % (type))
            if n > 1 and ((n % n_entries_int) == 1):
                fd.write("\n")
        fd.write(" }\n")

        atomic_positions_str = []
        for atom in structure:
            atomic_positions_str.append(
                "{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}\n".format(
                    coords=[atom.a, atom.b, atom.c]
                )
            )

        fd.write("xred {\n")
        fd.write("".join(atomic_positions_str))
        fd.write("}\n")

        fd.write(
            "acell {{ {acell[0]} {acell[0]} {acell[0]} }} \n".format(
                acell=[1 / Bohr]
            )
        )

        fd.write(
            "rprim {{ {cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n"
            "        {cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n"
            "        {cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}  }}\n"
            "".format(cell=structure.lattice.matrix)
        )

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

        # Obtain absorbing atom
        species = [structure[site].specie.symbol for site in sites]
        element = Element(species[0])

        cards = copy(self._cards)

        cards["edges"] = f"-{element.number} {self._edge}"
        # Estimate number of band
        if self._nbands_estimator == "heg":
            eRange = self._nbands_estimator_kwargs["eRange"]
            nbands = self._getCondBands(structure.lattice.volume, eRange)
            cards["nbands"] = -1 * nbands
        else:
            # if the _nbands_estimator != "heg", we expect the user to
            # probvide a number. Ocean accepts both positive and negative
            # values for nbands.
            try:
                nbands = int(self._nbands_estimator)
                cards["nbands"] = nbands
            except SyntaxError:
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
            cards["ngkpt"] = f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"
        else:
            raise ValueError("method for obtaining kmesh not supported")
        # Determine the diemac
        if self._bandgap is not None or self._diel is not None:
            if self._diel is not None:
                cards["diemac"] = self._diel
            elif self._bandgap is not None:
                cards["diemac"] = np.exp(3.5 / self._bandgap)
        elif self._bandgap is None and self._diel is None:
            if "bandgap" in kwargs.keys() and "diel" in kwargs.keys():
                bandgap = kwargs["bandgap"]
                diel = kwargs["diel"]
                if diel is not None:
                    if diel["poly_electronic"] is not None:
                        cards["diemac"] = diel["poly_electronic"]
                elif bandgap is not None:
                    if bandgap > 0.000001:
                        cards["diemac"] = np.exp(3.5 / bandgap)
                    else:
                        cards["diemac"] = 1000000
        # Determine the SCF? convergence threshold
        cards["toldfe"] = self._defaultConvPerAtom * len(structure)
        # Change screen.nkpt -Int to triplet of ints
        if "screen.nkpt" in cards and len(cards["screen.nkpt"].split()) == 1:
            cards["screen.nkpt"] = self._oceanKptSampling(
                structure.lattice.matrix, float(cards["screen.nkpt"])
            )

        # OCEAN will calculate every atom within the same specie
        path = target_directory  # / Path(f"{site:03}_{specie}")
        path.mkdir(exist_ok=True, parents=True)
        filepath_xas = path / "ocean.in"
        self._write_ocean_in(filepath_xas, structure, cards)

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

            with open(path / ("photon%d" % (photonCount)), "w") as f:
                f.write(mode + "\n")
                f.write("cartesian %f %f %f \n" % (dir1[0], dir1[1], dir1[2]))
                f.write("end\n")
                f.write("cartesian %f %f %f \n" % (dir2[0], dir2[1], dir2[2]))
                f.write("end\n")
                f.write("4966\n")
                # 4966 is hard coded for Ti
                # NEED TO FIX THIS (probably by moving it to a lookup table inside OCEAN)
                f.write(str(weight) + "\n")
                f.close

        return {"pass": True, "errors": dict()}

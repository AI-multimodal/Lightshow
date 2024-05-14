from pathlib import Path
from warnings import warn

from monty.json import MSONable
from pymatgen.core.periodic_table import Element

from lightshow.parameters._base import _BaseParameters

FDMNES_DEFAULT_CARDS = {
    "Energpho": True,
    "Memory_save": True,
    "Quadrupole": False,
    "Relativism": False,
    "Spinorbit": None,
    "SCF": True,
    "SCFexc": False,
    "Screening": False,
    "Full_atom": False,
    "TDDFT": False,
    "PBE96": False,
}


class FDMNESParameters(MSONable, _BaseParameters):
    """A one-stop-shop for all the different ways to modify input parameters
    for an FDMNES calculation.

    Parameters
    ----------
    cards : dict
        A dictionary of the cards to be control the parameters in the
        FDMNES calculations. The key of the dictionary corresponds to the
        parameters in FDMNES; the values are the correspongding values.
        In LightShow, recommended parameters to run Fdmnes is provided in
        ``FDMNES_DEFAULT_CARDS``, which looks something like

        .. code-block:: python

            cards = {
                "Energpho": True,
                "Memory_save": True,
                "Quadrupole": False,
                "Relativism": False,
                "Spinorbit": None,
                "SCF": True,
                "SCFexc": False,
                "Screening": False,
                "Full_atom": False,
                "TDDFT": False,
                "PBE96": False,
            }

        The detailed description of the FDMNES parameters can be find at
        its official website. If the user wants to change some parameters,
        they can just add a key-value pair to the cards.

    e_range : str
        The energy range E that one defines in the input is the energy of
        the photoelectron relative to the phonon level. If one wants the
        output energy relative to the Fermi level, put ``Energpho`` = False.
    edge : str
        The XAS edge of the calculation.
    radius : float
        FDMNES uses clusters for its calculations. The ``radius`` parameter
        determines how large to make the cluster. It is calculated from the
        absorbing atom center in units of Angstroms. The cluster radius is
        applicable to both SCF and XAS caluclations.
    name : str
        The name of the calculation.
    """

    def __init__(
        self,
        cards=FDMNES_DEFAULT_CARDS,
        e_range="-5. 0.5 60.",
        edge="K",
        radius=7.0,
        name=None,
    ):
        self._cards = cards
        self._radius = radius
        self._range = e_range

        self._name = name if name is not None else "FDMNES"
        self._edge = edge

        self.validate_edge()

    def validate_edge(self):
        """
        Validates and adjusts the edge attribute based on standard edge choices
        supported by FDMNES.

        Edge types recognized:
        - 'K', 'L1', 'L2', 'L3', 'L23'
        - 'M1', 'M2', 'M3', 'M23'
        - 'M4', 'M5', 'M45'
        - 'N1', 'N2', 'N3', 'N23'
        - 'N4', 'N5', 'N45'

        Warnings:
            Warns if the provided edge is 'L' that it is being modified to
            'L23'.
            Warns if the provided edge is not recognized and is being set to
            'K'.

        Returns:
            None.
        """

        valid_edges = [
            "K",
            "L1",
            "L2",
            "L3",
            "L23",
            "M1",
            "M2",
            "M3",
            "M23",
            "M4",
            "M5",
            "M45",
            "N1",
            "N2",
            "N3",
            "N23",
            "N4",
            "N5",
            "N45",
        ]
        if self._edge == "L":
            warn("Edge 'L' changed to 'L23' for FDMNES compatibility.")
            self._edge = "L23"
        elif self._edge not in valid_edges:
            warn(f"Edge {self._edge} not recognized. Defaulting to 'K'.")
            self._edge = "K"

    def get_FDMNESinput(self, structure, Z_absorber):
        """Constructs and returns a dictionary corresponds to the
        parameters in FDMNES optimized for the input structure and edge
        based on FDMNES documentation recommendations

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
            The Pymatgen structure. Note that the ``Z_absorber`` must
            correspond to the provided structure.
        Z_absorber : int
            Atomic number of the absorbing specie.

        Returns
        -------
        dict
            A dictionary of FDMNES input parameters.
        """

        cards = self._cards.copy()
        species_z_list = [species.Z for species in structure.species]

        transition_metal_ranges = [range(21, 31), range(39, 49), range(57, 81)]

        if self._edge == "K":
            if "Nonrelat" not in cards.keys():
                cards["Spinorbit"] = True
                warn(
                    "Spin-orbit has been turned on for K-edge calculation "
                    "for accuracy. The simulation is typically 4 to 8 times "
                    "longer and need 2 times more memory space. To turn"
                    " it off, set 'Nonrelat' = True. "
                )
            if any(Z_absorber in r for r in transition_metal_ranges):
                cards["Quadrupole"] = True

        elif self._edge == "L23" and Z_absorber in range(21, 26):
            cards["TDDFT"] = True

        if any(z > 36 for z in species_z_list):
            cards["Relativism"] = True

        if any(z > 50 for z in species_z_list):
            cards["Spinorbit"] = True

        if 8 in species_z_list:
            cards["Full_atom"] = True

        return cards

    def write(self, target_directory, **kwargs):
        """Writes the input files for the provided structure and absorber.
        In the case of Fdmnes, if Z_absorber is None, then the absorbing specie
        is the first one in the atom list.

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the FEFF input files.
        **kwargs
            Must contain the ``structure`` key (the
            :class:`pymatgen.core.structure.Structure` of interest) and the
            ``Z_absorber`` key (an int indicates the atomic number of the
            absorbing chemical specie).

        Returns
        -------
        dict
            A dictionary containing the status and errors key. In the case of
            FDMNES, there are no possible errors at this stage other than
            critical ones that would cause program termination, so the returned
            object is always
            ``{"pass": True, "errors": dict(), "path": ...}``.
        """

        structure = kwargs["structure"]
        Z_absorber = kwargs["Z_absorber"]

        if structure is None:
            raise ValueError("Structure must be provided.")

        # prepare lattice parameters, atomic numbers and fractional coordinates
        a, b, c = structure.lattice.abc
        alpha, beta, gamma = structure.lattice.angles
        atomic_numbers = structure.atomic_numbers
        scaled_positions = structure.frac_coords

        if Z_absorber is None:
            Z_absorber = atomic_numbers[0]
            warn(
                "Z_absorber is not provided, apply the first atom specie"
                f"to be the absorbing specie Z={atomic_numbers[0]} "
            )

        element_absorber = Element.from_Z(Z_absorber).symbol
        target_directory = Path(target_directory)
        target_directory.mkdir(exist_ok=True, parents=True)

        fdmnesinput = self.get_FDMNESinput(structure, Z_absorber)
        filepath = target_directory / f"{element_absorber}_in.txt"

        with open(filepath, "w") as f:
            f.write("Filout\n")
            f.write(f"  {element_absorber}\n\n")

            f.write("Range\n")
            f.write(f"  {self._range}\n\n")

            f.write("Radius\n")
            f.write(f"  {self._radius}\n\n")

            for key, value in fdmnesinput.items():
                if value:
                    f.write(f"{key}\n\n")

            f.write("Crystal \n")
            f.write(
                f"{a:.4f} {b:.4f} {c:.4f} {alpha:.1f} {beta:.1f} {gamma:.1f}\n"
            )
            for atomic_number, pos in zip(atomic_numbers, scaled_positions):
                f.write(
                    f"  {atomic_number} {pos[0]:.4f} {pos[1]:.4f} "
                    f"{pos[2]:.4f}\n"
                )

            f.write("\nZ_Absorber\n")
            f.write(f"  {Z_absorber}\n\n")
            f.write("Edge \n")
            f.write(f"  {self._edge}\n\n")
            f.write("Convolution \n\n")
            f.write("End")

        return {"pass": True, "errors": dict(), "path": str(filepath)}

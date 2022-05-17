from copy import copy
from functools import lru_cache

import numpy as np
from pathlib import Path
import warnings

from monty.json import MSONable
from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.inputs import Incar as pmgIncar
from pymatgen.io.vasp.inputs import Kpoints as pmgKpoints


class Incar(pmgIncar):
    """Inherits the pymatgen Incar class class. Allows for the standard
    manipulation of the INCAR file as a dictionary but with a few extra
    methods, such as saving and loading from file."""

    DEFAULT_SCF = {
        "ALGO": "Normal",
        "EDIFF": 1e-05,
        "IBRION": 2,
        "ISIF": 2,
        "ISMEAR": 0,
        "ISPIN": 2,
        "ISYM": 2,
        "KPAR": 1,
        "LAECHG": True,
        "LCHARG": True,
        "LORBIT": 11,
        "LREAL": "Auto",
        "LWAVE": False,
        "NCORE": 1,
        "NELM": 200,
        "NELMIN": 6,
        "NSIM": 16,
        "NSW": 0,
        "PREC": "Accurate",
        "SIGMA": 0.05,
        "SYMPREC": 1e-05,
    }

    DEFAULT_COREHOLE = {
        "ALGO": "Normal",
        "CH_LSPEC": True,
        "CH_NEDOS": 40000,
        "CH_SIGMA": 0.05,
        "CLL": 0,
        "CLN": 1,
        "CLNT": 1,
        "CLZ": 1.0,
        "EDIFF": 1e-05,
        "EDIFFG": -0.01,
        "ICORELEVEL": 2,
        "ISMEAR": 0,
        "ISPIN": 2,
        "KPAR": 1,
        "LAECHG": True,
        "LCHARG": False,
        "LREAL": "Auto",
        "LWAVE": False,
        "NCORE": 1,
        "NELM": 200,
        "NSIM": 16,
        "PREC": "Accurate",
        "SIGMA": 0.05,
        "SYMPREC": 1e-05,
    }

    def _adj_mag(
        self,
        struct,
        config_dict={
            "INCAR": {
                "LASPH": True,
                "MAGMOM": {
                    "Ce": 5,
                    "Ce3+": 1,
                    "Co": 0.6,
                    "Co3+": 0.6,
                    "Co4+": 1,
                    "Cr": 5,
                    "Dy3+": 5,
                    "Er3+": 3,
                    "Eu": 10,
                    "Eu2+": 7,
                    "Eu3+": 6,
                    "Fe": 5,
                    "Gd3+": 7,
                    "Ho3+": 4,
                    "La3+": 0.6,
                    "Lu3+": 0.6,
                    "Mn": 5,
                    "Mn3+": 4,
                    "Mn4+": 3,
                    "Mo": 5,
                    "Nd3+": 3,
                    "Ni": 5,
                    "Pm3+": 4,
                    "Pr3+": 2,
                    "Sm3+": 5,
                    "Tb3+": 6,
                    "Tm3+": 2,
                    "V": 5,
                    "W": 5,
                    "Yb3+": 1,
                },
            },
            "KPOINTS": {"reciprocal_density": 100},
        },
    ):
        """Adjusts the INCAR to use the default magmom and +U from pymatgen.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            The structure data for which the INCAR will be built
        config_dict : dict, optional
            The default settings for magmom from Materials Project
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            incar_tmp = DictSet(struct, config_dict).incar

        # Using spin unpolarized calculations for nm
        magmom_set = set(incar_tmp["MAGMOM"])
        if len(magmom_set) == 1 and 0 in magmom_set:
            incar_tmp["ISPIN"] = 1
            incar_tmp.pop("MAGMOM")

        # Assign the tags & values of magmoms and U to the INCAR
        for k, v in incar_tmp.items():
            self[k] = v

    def _adj_u(
        self,
        elements,
        config_dict={
            "Co": 3.32,
            "Cr": 3.7,
            "Fe": 5.3,
            "Mn": 3.9,
            "Mo": 4.38,
            "Ni": 6.2,
            "V": 3.25,
            "W": 6.2,
        },
    ):
        """Uses the Materials Project default values for the Hubburd U value,
        but since the LDAUU tag shoud match the elements in the POSCAR, weneed
        to add the u value (if needed) at the very end stage,  e.g. at the
        function write_single_VASP_files

        Parameters
        ----------
        elements : str
            string of element which should be the same as the element line
            in POSCAR
        config_dict : dict, optional
            default U values for the transition metal elements as in Materials Project
        """

        tm = ["Co", "Cr", "Fe", "Mn", "Mo", "Ni", "V", "W"]
        o_f = ["O", "F"]
        len_tm = len([ele for ele in tm if ele in elements])
        len_of = len([ele for ele in o_f if ele in elements])
        if len_tm > 0 and len_of > 0:
            self["LDAU"] = True
            self["LDAUJ"] = [0 for i in range(len(elements))]
            self["LDAUL"] = [2 for i in range(len(elements))]
            self["LDAUPRINT"] = 1
            self["LDAUTYPE"] = 2
            ldauu = []
            for ele in elements:
                if ele in config_dict:
                    ldauu.append(config_dict[ele])
                else:
                    ldauu.append(0)
            self["LDAUU"] = ldauu

    @classmethod
    def from_default(cls, scf=True):
        """A simple classmethod for loading the default specs.

        Parameters
        ----------
        scf : bool, optional
            Description

        Returns
        -------
        TYPE
            Description
        """

        if scf:
            klass = cls(Incar.DEFAULT_SCF)
        else:
            klass = cls(Incar.DEFAULT_COREHOLE)
        klass.check_params()
        return klass

    def write(self, filename):
        assert "INCAR" == Path(filename).name
        self.write_file(filename)


class Kpoints(pmgKpoints):
    @staticmethod
    def custom(supercell, cutoff=32.0, max_radii=50.0):
        """Customizes the KPOINTS file. For a kmesh sampling, e.g. [m, n, p],
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
        Kpoints
        """

        klist = dict()
        rlatt = np.array(supercell.lattice.reciprocal_lattice.abc)

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

        return Kpoints(kpts=(k,))

    def write_file(self, filename, assert_name=True):
        if assert_name:
            assert "KPOINTS" == Path(filename).name
        super().write_file(filename)


class PotcarGenerator(MSONable):
    """Summary"""

    # https://www.vasp.at/wiki/index.php/Available_PAW_potentials
    # using "Recommended potentials for GW/RPA calculations"

    DEFAULT_ELEMENT_MAPPING = {
        "Ti": "Ti_sv_GW",
        "Zn": "Zn_sv_GW",
        "O": "O_GW",
        "H": "H_GW",
        "He": "He_GW",
        "Li": "Li_sv_GW",
        "Be": "Be_sv_GW",
        "B": "B_GW",
        "C": "C_GW",
        "N": "N_GW",
        "F": "F_GW",
        "Ne": "Ne_GW",
        "Na": "Na_sv_GW",
        "Mg": "Mg_sv_GW",
        "Al": "Al_GW",
        "Si": "Si_GW",
        "P": "P_GW",
        "S": "S_GW",
        "Cl": "Cl_GW",
        "Ar": "Ar_GW",
        "K": "K_sv_GW",
        "Ca": "Ca_sv_GW",
        "Sc": "Sc_sv_GW",
        "V": "V_sv_GW",
        "Cr": "Cr_sv_GW",
        "Mn": "Mn_sv_GW",
        "Fe": "Fe_sv_GW",
        "Co": "Co_sv_GW",
        "Ni": "Ni_sv_GW",
        "Cu": "Cu_sv_GW",
        "Ga": "Ga_d_GW",
        "Ge": "Ge_d_GW",
        "As": "As_GW",
        "Se": "Se_GW",
        "Br": "Br_GW",
        "Kr": "Kr_GW",
        "Rb": "Rb_sv_GW",
        "Sr": "Sr_sv_GW",
        "Y": "Y_sv_GW",
        "Zr": "Zr_sv_GW",
        "Nb": "Nb_sv_GW",
        "Mo": "Mo_sv_GW",
        "Tc": "Tc_sv_GW",
        "Ru": "Ru_sv_GW",
        "Rh": "Rh_sv_GW",
        "Pd": "Pd_sv_GW",
        "Ag": "Ag_sv_GW",
        "Cd": "Cd_sv_GW",
        "In": "In_d_GW",
        "Sn": "Sn_d_GW",
        "Sb": "Sb_d_GW",
        "Te": "Te_GW",
        "I": "I_GW",
        "Xe": "Xe_GW",
        "Cs": "Cs_sv_GW",
        "Ba": "Ba_sv_GW",
        "La": "La_GW",
        "Ce": "Ce_GW",
        "Hf": "Hf_sv_GW",
        "Ta": "Ta_sv_GW",
        "W": "W_sv_GW",
        "Re": "Re_sv_GW",
        "Os": "Os_sv_GW",
        "Ir": "Ir_sv_GW",
        "Pt": "Pt_sv_GW",
        "Au": "Au_sv_GW",
        "Hg": "Hg_sv_GW",
        "Tl": "Tl_d_GW",
        "Pb": "Pb_d_GW",
        "Bi": "Bi_d_GW",
        "Po": "Po_d_GW",
        "At": "At_d_GW",
        "Rn": "Rn_d_GW",
    }

    def __init__(self, root, element_mapping=dict()):
        self._root = Path(root).resolve
        assert Path(self._root).exists()
        self._element_mapping = copy(
            PotcarGenerator.DEFAULT_ELEMENT_MAPPING
        ).update(element_mapping)

    def check_POTCAR_exists(self, list_of_elements):
        """Iterates through a list of element symbols and checks to see that a
        POTCAR file exists for each one of them.

        Parameters
        ----------
        list_of_elements : list
            A list of str, e.g., ``["Ti", "O"]``.

        Returns
        -------
        list
            A list of bool, where each entry corresponds to the true/false
            value of the POTCAR for that element existing.
        """

        el = []
        for element in list_of_elements:
            try:
                element_dir = self._element_mapping[element]
                path = Path(self._root) / Path(element_dir) / Path("POTCAR")
                el.append(path.exists())
            except KeyError:
                el.append(False)
        return el

    @lru_cache(256)
    def _get_element_lines(self, element):

        # Get the element directory using the element as a key
        assert isinstance(element, str)
        element_dir = self._element_mapping[element]

        # Construct the path
        path = Path(self._root) / Path(element_dir) / Path("POTCAR")

        # Load in the file
        with open(path, "r") as f:
            lines = f.readlines()

        return lines

    @lru_cache(256)
    def get_n_valence_electrons(self, element):
        """Reads the total number of valence electrons from the second line of
        the potential file.

        Parameters
        ----------
        element : str
            The element, e.g., ``"Ti"``.

        Returns
        -------
        int
            The number of valence electrons for that atom as specified by the
            potential file.
        """

        element_dir = self._element_mapping[element]

        # Construct the path
        path = Path(self._root) / Path(element_dir) / Path("POTCAR")

        # Load in the file
        with open(path, "r") as f:
            _ = f.readline()
            nelec = f.readline()

        return int(float(nelec.strip()))

    def get_total_valence_electrons(self, structure):
        """Gets the total number of valence electrons in a structure by
        counting the number of valence electrons per element.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure

        Returns
        -------
        int
            The total valence electrons in that structure. Note it counts the
            total valence in the full structure provided, so supercells will
            have more total valence electrons than primitive cells.
        """

        elements = [str(s.specie) for s in structure]
        return int(sum([self.get_n_valence_electrons(ii) for ii in elements]))

    # TODO: Fanchen: can you fix the equation in this docstring? Doesn't look
    # correct, or at least it's not clear
    def get_total_bands_via_heg_approximation(self, structure, eRange):
        """Gets an approximation for the total number of bands to use in the
        calculation based on the free electron gas model, such that a range of
        ``eRange`` is covered in the spectrum. The method assumes all electrons
        are free electrons, including those in the valence bands. The total
        number of electrons will then be calculated as
        :math:`N = N_\\mathrm{v} + N_\\mathrm{c}`, where :math:`N_\\mathrm{v}`
        is the total number of valence electrons and
        :math:`N_\\mathrm{c}` is the number of conduction band electrons above
        the Fermi energy as specified by the provided energy range. The final
        expression (in atomic units) is

        .. math::

            E_\\mathrm{range}/c_1 = c_2 / (V/N_vc * 3/ (4 pi))^2/3
            - c_2 / (V/N_v * 3 / (4 pi))^2/3 N_v

        where :math:`c_1 = 27.2114` and :math:`c_2 = 1.841`,
        is summed over all valence electrons obtained from POTCAR.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure
        eRange : float
            The approximate energy range desired in the spectrum, starting at
            the Fermi energy.

        Returns
        -------
        int
            Approximation for the number of bands to use in computing the XANES
            spectrum as computed using the (3D) homogeneous electron gas
            approximation.
        """

        ANG2BOHR = 1.88873
        volume_bohr = structure.volume * ANG2BOHR**3
        Nv = self.get_total_valence(structure)

        tmp1 = eRange / 27.2114 + 1.841 / (
            volume_bohr / Nv * 3.0 / 4.0 / np.pi
        ) ** (2 / 3)
        tmp2 = (1.841 / tmp1) ** (2 / 3) / 3.0 * 4.0 * np.pi
        Nvc = volume_bohr / tmp2 * 0.65
        return int(np.ceil(Nvc))

    def write(self, target_directory, elements):
        """Writes the POTCAR file to the path provided by stacking the
        potentials for the provided elements.

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the POTCAR file.
        elements : list
            A list of str: the element symbols to use.
        """

        path = Path(target_directory) / Path("POTCAR")

        lines = []
        for element in elements:
            lines.extend(self._get_element_lines(element))

        with open(path, "w") as f:
            for line in lines:
                f.write(line)

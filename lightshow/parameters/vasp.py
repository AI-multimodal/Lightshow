from copy import copy, deepcopy
from functools import lru_cache
from pathlib import Path
import warnings

from monty.json import MSONable
import numpy as np
from pymatgen.io.vasp.sets import DictSet
from pymatgen.io.vasp.inputs import Incar as pmgIncar
from pymatgen.io.vasp.inputs import Kpoints as pmgKpoints
from pymatgen.io.vasp.inputs import Poscar as pmgPoscar

from lightshow.parameters._base import _BaseParameters


class Incar(pmgIncar):
    """Inherits the Pymatgen Incar class. Allows for the standard manipulation
    of the INCAR file as a dictionary but with a few extra methods."""

    DEFAULT_NEUTRAL = {
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

    def adj_mag(
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
        """Adjusts the INCAR to use the default magmom and +U from Pymatgen.

        Parameters
        ----------
        struct : pymatgen.core.structure.Structure
            The structure data for which the INCAR will be built.
        config_dict : dict, optional
            The default MAGMOM settings from the Materials Project.
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

    def adj_u(
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
        """Updates the INCAR using the Materials Project default values for the
        Hubburd U values.

        Parameters
        ----------
        elements : list
            A list of str, containing the elements in the same order as the
            element line in the POSCAR.
        config_dict : dict, optional
            Default Materials Project U values for the transition metal
            elements.
        """

        tm = list(config_dict.keys())
        o_f = ["O", "F"]  # TODO: Fanchen, why just O and F?
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
    def from_default(cls, neutral=True):
        """A simple classmethod for loading the default specs.

        Parameters
        ----------
        neutral : bool, optional
            If True, returns the default for a neutral potential electronic
            structure relaxation. Otherwise, assumes a corehole potential.

        Returns
        -------
        Incar
        """

        if neutral:
            klass = cls(Incar.DEFAULT_NEUTRAL)
        else:
            klass = cls(Incar.DEFAULT_COREHOLE)
        return klass

    def check_params(self):
        """Raises warnings for erroneous INCAR tags/parameters (see the
        `pymatgen.io.vasp.inputs documentation <https://pymatgen.org/
        pymatgen.io.vasp.inputs.html>`_ for precisely what this does). In
        addition, will actually raise a ValueError if "NBANDS" is not present
        in the INCAR.

        Raises
        ------
        ValueError
            If "NBANDS" is not in the INCAR ``keys``.
        """

        if "NBANDS" not in self.keys():
            raise ValueError("NBANDS not present in INCAR parameters")
        super().check_params()

    def write(self, target_directory):
        """Writes the INCAR file to the target directory.

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the INCAR file.
        """

        self.write_file(Path(target_directory) / Path("INCAR"))


class Kpoints(pmgKpoints):
    """Method for constructing the Kpoints objects.

    .. note::

        There are many more custom ways to instantiate the Kpoints objects than
        just the :class:`Kpoints.custom` method listed below. The full list
        can be found in the `Pymatgen documentation <https://pymatgen.org/
        pymatgen.io.vasp.inputs.html#pymatgen.io.vasp.inputs.Kpoints>`_.
    """

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

    def write(self, target_directory):
        """Writes the KPOINTS file to the target directory.

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the KPOINTS file.
        """

        self.write_file(Path(target_directory) / Path("KPOINTS"))


class PotcarConstructor(MSONable):
    """A helper class used to create the POTCAR input files for VASP
    calculations. The default choices for the types of potentials chosen is
    detailed in :class:`.PotcarConstructor.DEFAULT_ELEMENT_MAPPING` and is
    chosen from the recommended potentials for GW/RPA calculations found
    `here <https://www.vasp.at/wiki/index.php/
    Available_PAW_potentials#Recommended_potentials_for_GW.
    2FRPA_calculations>`_.

    .. warning::

        The potential files fall under the VASP license. If you do not have a
        VASP license, you cannot use VASP potential files. Please take a look
        at the VASP `website <https://www.vasp.at/>`_ for more details.

    Parameters
    ----------
    root : str
        The location of the potential file-containing directories, e.g.,
        ``"Ti_sv_GW"``.
    element_mapping : dict
        Overrides specific defaults in the
        :class:`.PotcarConstructor.DEFAULT_ELEMENT_MAPPING`. Note that only
        provided values will be overridden, and others will fallback to the
        defaults.
    """

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
        self._root = Path(root).resolve()
        assert Path(self._root).exists()
        self._element_mapping = copy(PotcarConstructor.DEFAULT_ELEMENT_MAPPING)
        self._element_mapping.update(element_mapping)

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

    def check_POSCAR_valid(self, poscar):
        """Checks the poscar against the potcar data to ensure that all elements
        have a valid potential.

        Parameters
        ----------
        poscar : Poscar

        Returns
        -------
        str
            The elements that do not exist, else None if no error.
        """

        s = poscar.site_symbols
        elements_exist = self.check_POTCAR_exists(s)

        if not all(elements_exist):
            dont_exist = ",".join(
                [
                    xx
                    for xx, tf in zip(poscar.site_symbols, elements_exist)
                    if not tf
                ]
            )
            return dont_exist

        return None

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
        """For a given element, reads the total number of valence electrons
        from the second line of the potential file.

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
    def get_total_bands_via_heg(self, structure, eRange):
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
        Nv = self.get_total_valence_electrons(structure)

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
            A list of str: the elements to use in the POTCAR file, in the
            correct order.
        """

        path = Path(target_directory) / Path("POTCAR")

        lines = []
        for element in elements:
            lines.extend(self._get_element_lines(element))

        with open(path, "w") as f:
            for line in lines:
                f.write(line)


class Poscar(pmgPoscar):
    """Custom Poscar file, lightweight wrapper for the
    :class:`pymatgen.io.vasp.inputs.Poscar` class.

    .. warning::

        This :class:`.Poscar` class is only tested with VASP >6.2. It is likely
        compatible with earlier versions, but the user should proceed with
        caution in those cases.
    """

    def __len__(self):
        return len(self.structure)

    def write_single_site(
        self, target_directory, site_index=None, check_atom_type=None
    ):
        """Writes the VASP POSCAR file for the specified ``site_index``.

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the POSCAR file.
        site_index : int, optional
            In corehole calculations, we must move the site of interest to the
            first entry in the POSCAR file, and ensure that the element lines
            in said POSCAR are changed accordingly. For example,
            ``Ti O -> Ti Ti O`` and ``6 10 -> 1 5 10``, are changed in order to
            properly reflect the new ordering of atoms. This is to simplify how
            the POTCAR file is saved, and indicates for VASP to use a different
            potential for the atom with the corehole on it. Default is None
            (indicating no changes to be made, this is useful for neutral
            potential calculations).
        check_atom_type : str, optional
            If not None, asserts that the found atom type in the POSCAR lines
            is the same as the provided atom type.

        Returns
        -------
        list
            A list of str: the ordering of the elements for the POTCAR file.
        """

        path = Path(target_directory) / Path("POSCAR")

        # This might be specific to certain versions of VASP
        STRUCTURE_START_INDEX = 8

        # In the case of e.g. the neutral potential calculation, the base class
        # method write_file can be used
        if site_index is None:
            self.write_file(path)
            return self.site_symbols

        # Get the default lines for the POSCAR file
        lines = self.get_string().split("\n")

        # We need to make some modifications depending on the type of
        # calculation this is. Lines indexed by 5 and 6 are the ones that
        # need to be changed, and lines 8 onward contain the structure data
        # to be swapped.
        assert isinstance(site_index, int)
        assert site_index >= 0
        site_index += STRUCTURE_START_INDEX

        # Get the atom type of the site provided
        atom_type = lines[site_index].split()[-1]

        # Assert if provided
        if check_atom_type is not None:
            assert atom_type == check_atom_type

        # Get the location of the atom type
        atom_type_loc = np.where(np.array(self.site_symbols) == atom_type)[
            0
        ].item()
        new_line_5 = [atom_type, *self.site_symbols]
        line_6_asint = [int(xx) for xx in self.natoms]
        line_6_asint[atom_type_loc] -= 1
        new_line_6_asint = [1, *line_6_asint]
        new_line_6 = [str(xx) for xx in new_line_6_asint]

        # Use the new lines 5 and 6
        lines[5] = " ".join(new_line_5)
        lines[6] = " ".join(new_line_6)

        # Execute the swap. First, get the line of interest
        move_to_front = lines[site_index]

        # Delete that entry from the list
        del lines[site_index]

        # And insert back at the top
        lines.insert(STRUCTURE_START_INDEX, move_to_front)

        # And save
        with open(path, "w") as f:
            for line in lines:
                f.write(f"{line}\n")

        return new_line_5


class VASPParameters(MSONable, _BaseParameters):
    """A one-stop-shop for modifying the VASP parameters for some calculation.
    Due to the relative complexity of setting up a VASP calculation (as opposed
    to e.g. FEFF), the :class:`.VASPParameters` object...

    Parameters
    ----------
    incar : dict or :class:`.Incar`
        The INCAR VASP parameter information. If of type dict, then an
        :class:`.Incar` class will be instantiated from it. Can also be of type
        :class:`.Incar` directly. The number of bands does not have to be
        directly specified, as it can be estimated based on the structural
        information as long as ``nbands_estimator`` is not ``None``.
    kpoints_method : str
        Precise name of the classmethod defined in
        :class:`pymatgen.io.vasp.inputs.Kpoints`, or ``"custom"``, which is
        defined in Lightshow's :class:`.Kpoints` module.
    potcar_directory : os.PathLike
        The location in which the potential files are stored. These files
        should be stored in a precise directory format, specifically: the
        first level in should contain potential files of the form
        ``Pb_d_GW``, ``Pb_sv_GW``, ``Pd``, etc. For some element, the specific
        choice is determined first by the ``potcar_element_mapping`` and then
        by the defaults in the :class:`.PotcarConstructor` class.
    kpoints_method_kwargs : dict
        Optional arguments to pass to the classmethod used to construct the
        kpoints.
    potcar_element_mapping : dict
        A custom element mapping which will override, element-by-element,
        defaults in the :class:`.PotcarConstructor` class.
    nbands_estimator : str
    """

    @property
    def incar(self):
        return self._incar

    def get_KPOINTS(self, structure):
        """Calls the static methods defined on Kpoints (and it's corresponding
        base class) using the set arguments.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure

        Returns
        -------
        pymatgen.io.vasp.inputs.Kpoints
        """

        f = getattr(Kpoints, self._kpoints_method)
        return f(structure, **self._kpoints_method_kwargs)

    def estimate_n_bands(self, structure):
        """Uses the ``nbands_estimator`` information provided at instantiation
        to estimate the number of bands to use for a calculation.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure

        Returns
        -------
        int

        Raises
        ------
        KeyError
            If ``nbands_estimator`` is None, and there is no ``NBANDS`` key
            provided in the INCAR. Will also throw a KeyError if the correct
            parameter information (for the specified method) is not contained
            in ``nbands_parameters``.
        ValueError
            If the specified ``nbands_method`` is not one of ``None``,
            ``int`` or ``"heg"``.
        """

        cond2 = isinstance(self._nbands_estimator, int)
        if self._nbands_estimator is None or cond2:
            if self._incar.get("NBANDS") is None:
                raise KeyError(
                    "nbands_estimator is None or int and the INCAR does not "
                    "contain  a specified number of bands"
                )

        if self._nbands_estimator is None:
            return self._incar["NBANDS"]

        if cond2:
            nvalence = self._potcar_constructor.get_n_valence_electrons(
                structure
            )
            nb = nvalence * self._nbands_estimator
        elif self._nbands_estimator == "heg":
            nb = self._potcar_constructor.get_total_bands_via_heg(
                structure, self._nbands_parameters["erange"]
            )
        else:
            raise ValueError(f"Unknown nbands method={self._nbands_estimator}")

        return nb

    @property
    def name(self):
        return self._name

    def __init__(
        self,
        incar,
        potcar_directory,
        kpoints_method="custom",
        kpoints_method_kwargs={"cutoff": 32.0, "max_radii": 50.0},
        potcar_element_mapping=dict(),
        nbands_estimator=None,
        nbands_parameters=dict(),
        max_bands=int(1e16),
        force_spin_unpolarized=False,
        name="VASP",
    ):
        # Load the INCAR information
        if isinstance(incar, Incar):
            self._incar = incar
        elif isinstance(incar, dict):
            self._incar = Incar(incar)
        else:
            raise ValueError(
                f"Unknown incar type: {type(incar)}. Should be of type "
                "Incar or dict"
            )

        # Method for generating the KPOINTS file
        self._kpoints_method = kpoints_method
        self._kpoints_method_kwargs = kpoints_method_kwargs

        # POTCAR information
        self._potcar_directory = potcar_directory
        self._potcar_element_mapping = potcar_element_mapping
        self._potcar_constructor = PotcarConstructor(
            potcar_directory, potcar_element_mapping
        )

        # Estimator method for the number of bands
        self._nbands_estimator = nbands_estimator
        self._nbands_parameters = nbands_parameters

        # Other parameters
        self._max_bands = max_bands
        self._force_spin_unpolarized = force_spin_unpolarized
        self._name = name

    def write(self, target_directory, **kwargs):
        """Summary

        Parameters
        ----------
        target_directory : os.PathLike
            The target directory to which to save the VASP input files.
        **kwargs
            Must contain the ``structure_sc`` key (the
            :class:`pymatgen.core.structure.Structure` of interest) and the
            ``sites`` key (a list of int, where each int corresponds to the
            site index of the site to write).

        Returns
        -------
        dict
            A dictionary containing the status and errors key.

        Raises
        ------
        ValueError
            If the symmetrically inequivalent sites provided is None but the
            CH_LSPEC flag is True in the INCAR.
        """

        errors = dict()

        incar = deepcopy(self._incar)
        sites = kwargs.get("sites")
        ch_lspec = incar.get("CH_LSPEC", False)
        if sites is None and ch_lspec:
            raise ValueError(
                "Corehole calculation has been selected via the presence of "
                "CH_LSPEC in the INCAR, but no sites were provided."
            )

        structure = deepcopy(kwargs["structure_sc"])

        # Check the validity of the calculations
        nb = self.estimate_n_bands(structure)
        if nb > self._max_bands:
            errors["n_bands"] = nb

        poscar = Poscar(structure)
        potcar_check = self._potcar_constructor.check_POSCAR_valid(poscar)
        if potcar_check is not None:
            errors["potentials_do_not_exist"] = potcar_check

        incar["NBANDS"] = self.estimate_n_bands(structure)
        if self._force_spin_unpolarized:
            incar["ISPIN"] = 1
            try:
                incar.pop("MAGMOM")
            except KeyError:
                pass
        else:
            incar.adj_mag(structure)
        incar.check_params()

        if len(errors) > 0:
            return {"pass": False, "errors": errors}

        kpoints = self.get_KPOINTS(structure)
        species = [structure[site].specie.symbol for site in sites]

        # Green light for corehole calculations
        if ch_lspec:
            for site, specie in zip(sites, species):
                path = target_directory / Path(f"{site:03}_{specie}")
                path.mkdir(exist_ok=True, parents=True)
                elements = poscar.write_single_site(
                    path, site_index=site, check_atom_type=specie
                )
                self._potcar_constructor.write(path, elements)
                kpoints.write(path)
                incar.adj_u(elements)
                incar.write(path)

        # Else assume this is just a standard neutral potential calculation
        else:
            path = target_directory / Path("SCF")
            path.mkdir(exist_ok=True, parents=True)
            elements = poscar.write_single_site(
                path, site_index=None, check_atom_type=None
            )
            self._potcar_constructor.write(path, elements)
            kpoints.write(path)
            incar.adj_u(elements)
            incar.write(path)

        return {"pass": True, "errors": dict()}

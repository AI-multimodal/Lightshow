"""Contains all relevant classes and functions for constructing databases of
various materials and clusters. This includes pulling and processing data from
the `Materials Project Database <https://materialsproject.org/>`_, as well as
utilizing existing data the user may have on their hard drive."""

import json
from datetime import datetime
from pathlib import Path
from shutil import copy2
from warnings import warn

from monty.json import MSONable
from mp_api.client import MPRester
from pymatgen.core.structure import Molecule, Structure
from tqdm import tqdm

from lightshow import pymatgen_utils
from lightshow.utils.environ_utils import get_API_key_from_environ


def _get_api_key(api_key):
    if api_key is None:
        api_key = get_API_key_from_environ()
    if api_key is None:
        raise ValueError(f"Invalid API key {api_key}")
    return api_key


def _get_method(method, mpr):
    """Get all the available search methods; if the given method is not
    present, the search will be performed using mpr.materials.search"""

    methods = [
        methods
        for methods in dir(mpr.materials)
        if methods is not callable and not methods.startswith("_")
    ]

    if method is None:
        print("Searching Materials Project with default method")
    elif method not in methods:
        warn(
            f"Provided method={method} not in available methods {methods}",
            "falling back on default search method.",
        )
        method = None

    return method


class Database(MSONable):
    """Contains all materials and metadata for some database."""

    @classmethod
    def from_files_molecule(
        cls,
        root,
        filename="*.xyz",
        lattice=None,
        pbar=True,
    ):
        """Searches for files matching the provided ``filename``, and assumes
        those files are structural files in a format compatible with
        ``Molecule.from_file``.

        Parameters
        ----------
        root : str
            The directory in which to begin the search.
        filename : str, optional
            The files to search for. Uses ``rglob`` to recursively find any
            files matching ``filename`` within the provided directory.
        lattice : list of floats, optional
            Lattice parameter used to construct the crystal lattice. If not
            provided, defaults to [20.0, 20.0, 20.0] Angstroms.
        pbar : bool, optional
            If True, will show a tqdm progress bar.

        Returns
        -------
        Database
        """

        if lattice is None:
            lattice = [20.0, 20.0, 20.0]

        structures = {}
        metadata = {}
        for key, path in enumerate(
            tqdm(Path(root).rglob(filename), disable=not pbar)
        ):
            key = f"{key:08}"
            molecule = Molecule.from_file(path)
            structures[key] = molecule.get_boxed_structure(*lattice)
            metadata[key] = {"origin": str(path)}
        return cls(structures=structures, metadata=metadata, supercells=dict())

    @classmethod
    def from_files(cls, root, filename="CONTCAR", pbar=True):
        """Searches for files matching the provided ``filename``, which can
        include wildcards, and assumes those files are structural files in a
        format that can be processed by ``Structure.from_file``. Each structure
        is given its own index, with the origin path stored in its metadata.

        .. code-block:: python

            {
                "0": struct1,
                "1": struct2,
                ...
            }

        Parameters
        ----------
        root : str
            The directory in which to begin the search.
        filename : str, optional
            The files to search for. Uses ``rglob`` to recursively find any
            files matching ``filename`` within the provided directory.
        pbar : bool, optional
            If True, will show a tqdm progress bar.

        Returns
        -------
        Database
        """

        structures = {}
        metadata = {}
        for key, path in enumerate(
            tqdm(Path(root).rglob(filename), disable=not pbar)
        ):
            key = f"{key:08}"
            struct = Structure.from_file(path)
            structures[key] = struct.get_primitive_structure()
            metadata[key] = {"origin": str(path)}
        return cls(structures=structures, metadata=metadata, supercells=dict())

    @classmethod
    def from_materials_project(cls, api_key=None, method=None, **kwargs):
        """Constructs the :class:`.Database` object by pulling structures and
        metadata directly from the Materials Project. This is a simple
        passthrough method which utilizes the MPRester.materials.search
        API of the Materials Project v2 API.

        Parameters
        ----------
        api_key : None, optional
            API key which can either be provided directly or is read from
            the MP_API_KEY environment variable.
        method : None, optional, str
            Keyword to get different information about materials'
            for e.g. 'thermo', 'xas', 'summary' etc. fetch information on
            thermodynamic properties, computed XAS data, large amount of
            amalgated data about the material, respectively.
            See https://api.materialsproject.org/docs for more details.
        **kwargs
            Description

        Returns
        -------
        Database
        """

        api_key = _get_api_key(api_key)
        method = kwargs.get("method")

        with MPRester(api_key) as mpr:
            method = _get_method(method, mpr=mpr)
            if method is not None:
                searched = getattr(mpr.materials, method).search(**kwargs)
            else:
                searched = mpr.materials.search(**kwargs)

        structures = {
            s.material_id.string: s.structure
            if hasattr(s, "structure")
            else None
            for s in searched
        }
        metadata = {s.material_id.string: s.dict() for s in searched}

        return cls(structures=structures, metadata=metadata, supercells=dict())

    def initialize_supercells(self, supercell_cutoff=9.0):
        """Initializes the supercells from the structures pulled from the
        Materials project.

        Parameters
        ----------
        supercell_cutoff : float, optional
            Parameter used for constructing the supercells. Default is 9
            Angstroms. TODO: this needs clearer documentation.
        """

        for key, prim in self._structures.items():
            supercell = pymatgen_utils.make_supercell(prim, supercell_cutoff)
            self._supercells[key] = supercell.get_sorted_structure()
        self._supercell_cutoff = supercell_cutoff
        self._supercells_initialized = True

    def initialize_inequivalent_sites(self):
        """Iterates through the structures and updates the metadata with
        keys corresponding to the inequivalent sites in the structure. This
        also tracks the atom types corresponding to the inequivalent sites and
        their multiplicities in the structure."""

        for key, prim in self._structures.items():
            info = pymatgen_utils.get_inequivalent_site_info(prim)
            self._metadata[key]["primitive"] = info

        # Do the same for the supercells
        for key, supercell in self._supercells.items():
            info = pymatgen_utils.get_inequivalent_site_info(supercell)
            mapping = pymatgen_utils.get_supercell_indexes_matching_primitive(
                self._structures[key], supercell, compare=30, r=20.0
            )
            info["prim-supercell-mapping"] = mapping
            self._metadata[key]["supercell"] = info

        self._inequivalent_sites_initialized = True

    @property
    def structures(self):
        """A dictionary of :class:`pymatgen.core.structure.Structure` objects.
        Contains the primitive structures. The keys are the IDs of the
        structure. In the case of data pulled from the Materials Project, these
        are the MPIDs. Otherwise, they are simply strings encoding some
        information about the origins of the structures.

        Returns
        -------
        dict
        """

        return self._structures

    @property
    def supercells(self):
        """A dictionary of :class:`pymatgen.core.structure.Structure` objects
        containing supercells and the same keys as ``structures``.

        Returns
        -------
        dict
        """

        if not self._supercells_initialized:
            warn("Run initialize_supercells(); supercells is probably empty")

        return self._supercells

    @property
    def metadata(self):
        """A dictionary of metadata information about the structure. This
        information is filled by the Materials Project when data is pulled,
        and might be empty if there is no associated metadata (in cases of
        using data already on disk, for example).

        Returns
        -------
        dict
        """

        return self._metadata

    @property
    def database_status(self):
        """A dictionary containing the current status of the database.
        Basically everything except the structures, metadata and supercells.

        Returns
        -------
        dict
        """

        return {
            key: value
            for key, value in self.as_dict().items()
            if key not in ["structures", "metadata", "supercells"]
        }

    def __init__(
        self,
        structures,
        metadata=None,
        supercells=None,
        supercell_cutoff=None,
        inequivalent_sites_initialized=False,
        supercells_initialized=False,
    ):
        """Initializer for the Database class. Note it is recommended to use
        the classmethods to initialize this object."""

        self._structures = structures
        if metadata is None:
            self._metadata = {}
        else:
            self._metadata = metadata
        if supercells is None:
            self._supercells = {}
        else:
            self._supercells = supercells
        self._supercell_cutoff = supercell_cutoff
        self._inequivalent_sites_initialized = inequivalent_sites_initialized
        self._supercells_initialized = supercells_initialized

    def _setup_preliminary_attributes(self):
        """Initializes supercells and inequivalent site info in the metadata
        if they're not already."""

        if not self._supercells_initialized:
            warn(
                "Initializing supercells with supercell_cutoff=9.0. "
                "Run initialize_supercells(supercell_cutoff=...) to set "
                "this cutoff manually."
            )
            self.initialize_supercells(9.0)

        if not self._inequivalent_sites_initialized:
            self.initialize_inequivalent_sites()

    @staticmethod
    def _get_site_indexes_matching_atom(info, species):
        """For an info dictionary of the form
        {
            "sites": inequivalent_sites,
            "species": species,
            "multiplicities": multiplicities
        }
        gets the inequivalent indexes matching the provided atom type/species.
        Also returns None if ``species`` is None.

        Parameters
        ----------
        info : dict
        species : str
            The atom type, e.g. "Ti", "O", etc.

        Returns
        -------
        list
        """

        if species is None:
            return None

        return [
            index
            for index, specie in zip(info["sites"], info["species"])
            if specie == species
        ]

    def _write_unit_cells(self, root, pbar=False):
        """A helper method for writing the unit cells in POSCAR format. This
        is convenient as the atom indexes saved as Lightshow runs correspond
        to this unit cell.

        Parameters
        ----------
        root : os.PathLike
        pbar : bool, optional
        """

        for key, structure in tqdm(self._structures.items(), disable=not pbar):
            fname = Path(root) / key / "POSCAR"
            structure.to(fmt="POSCAR", filename=str(fname))

    def _write_multiplicity(self, root, pbar=False):
        """Helper method for writing the site multiplicity information of
        the structure."""

        for key, metadata in tqdm(self._metadata.items(), disable=not pbar):
            fname = Path(root) / key / "multiplicity.json"
            multiplicities = {
                i_site: n_multi
                for i_site, n_multi in zip(
                    metadata["primitive"]["sites"],
                    metadata["primitive"]["multiplicities"],
                )
            }
            with open(fname, "w") as outfile:
                json.dump(multiplicities, outfile, indent=4, sort_keys=True)

    def write(
        self,
        root,
        absorbing_atoms=None,
        options=[],
        pbar=True,
        copy_script=None,
        write_unit_cells=True,
        write_multiplicity=True,
    ):
        """The core method of the :class:`.Database` class. This method will
        write all input files specified in the ``options`` parameter to disk.
        Of particular note is the directory structure, which is always
        consistent regardless of the type of calculation: At the first level is
        the ``key`` (usually an mpid) indexing the material of interest. At the
        next level is the user-specified "name" of the calculation, which are
        the first elements of each tuple in the ``options``. Next are
        the material-specific calculations. For example, in spectroscopy
        calculations, each absorbing atom will have it's own directory at this
        level. Within each of those are the input files for that calculation.
        For example::

            mp-390
                VASP
                    000_Ti
                    SCF
                FEFF-XANES
                    000_Ti
                        # ... (input files)
            mvc-11115
                # ...

        .. note::

            At the end of writing the files, a ``writer_metadata.json`` file
            will be saved along with the directories containing the input
            files. This metadata file contains all information about the
            parameters used to construct the input files, including those
            passed as arguments to this method.

        Parameters
        ----------
        root : os.PathLike
            The target root directory to save all of the input files.
        absorbing_atoms : str or list, optional
            The absorbing atom type symbol(s), e.g. ``"Ti"``. Note that if None,
            any calculations in which the absorbing atom is required (e.g. all
            spectroscopy) will be skipped. Only calculations that do not
            require absorbing atoms to be specified (e.g. neutral potential
            VASP electronic structure self-consistent procedure) will be
            performed. Note this can also be ``"all"``, in which case, every
            atom in the structure will have input files written for it.
        options : list, optional
            A list of :class:`lightshow.parameters._base._BaseParameters`
            objects or derived instances. The choice of options not only
            specifies which calculations to setup, each of the options also
            contains the complete set of parameters necessary to characterize
            each individual set of input files (for e.g. FEFF, VASP, etc.).
        pbar : bool, optional
            If True, enables the :class:`tqdm` progress bar.
        copy_script : os.PathLike
            If not ``None``, will copy the script in the provided path to each
            of the input file locations.
        write_unit_cells : bool, optional
            If True, writes the unit cells in the materials directory in
            POSCAR format. Very useful!
        write_multiplicity : bool, optional
            If True, writes the multiplicities of each atom in the unit cell
            to a multiplicities.json file.
        """

        self._setup_preliminary_attributes()

        root = str(Path(root).resolve())  # Get absolute path

        _write_all_atoms = False
        if absorbing_atoms == "all":
            _write_all_atoms = True
        elif not isinstance(absorbing_atoms, list):
            absorbing_atoms = [absorbing_atoms]

        writer_metadata = {
            "locals": {
                key: value
                for key, value in locals().items()
                if key not in ["self", "options"]
            },
            "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "options": [option.as_dict() for option in options],
            **self.database_status,
        }

        root = Path(root)
        root.mkdir(exist_ok=True, parents=True)
        writer_metadata_path = root / Path("writer_metadata.json")

        for key, supercell in tqdm(self._supercells.items(), disable=not pbar):
            primitive_cell = self._structures[key]
            primitive_info = self._metadata[key]["primitive"]
            supercell_info = self._metadata[key]["supercell"]
            prim_to_sc_index_mapping = supercell_info["prim-supercell-mapping"]

            if _write_all_atoms:
                absorbing_atoms = list(
                    set([s.specie.symbol for s in supercell])
                )

            for absorbing_atom in absorbing_atoms:
                # Gracefully skip any atoms that are not present in the
                # structure
                if not pymatgen_utils.atom_in_structure(
                    absorbing_atom, supercell
                ):
                    continue

                # If inequiv is None, that means that the absorbing_atom was not
                # specified (absorbing_atom is None)
                inequiv = self._get_site_indexes_matching_atom(
                    primitive_info, absorbing_atom
                )

                # If the for VASP and XSpectra calculations, use supercell;
                # otherwise, use unit cell structure
                # test if band_gap and diel in the self._metadata[key].keys()
                # if yes, read the bandgap and diel for OCEAN
                # if no, ignore them
                kwargs = {
                    "structure_sc": supercell,
                    "structure_uc": primitive_cell,
                    "sites": inequiv,
                    "index_mapping": prim_to_sc_index_mapping,
                }
                if key in self._metadata.keys():
                    if (
                        "band_gap" in self._metadata[key].keys()
                        and "diel" in self._metadata[key].keys()
                    ):
                        kwargs["bandgap"] = self._metadata[key]["band_gap"]
                        kwargs["diel"] = self._metadata[key]["diel"]

                # Write the files that we can
                for option in options:
                    path = root / Path(key) / Path(option.name)
                    status = option.write(path, **kwargs)
                    if not status["pass"]:
                        d = {"name": key, **status["errors"]}
                        warn(f"error: {key}+{option.name}: {d}")
                    if copy_script is not None and "paths" in status.keys():
                        for p in status["paths"]:
                            copy2(copy_script, p)

        if write_unit_cells:
            self._write_unit_cells(root, pbar=pbar)
        if write_multiplicity:
            self._write_multiplicity(root, pbar=pbar)

        # Save a metadata file (not a serialized version of this class) to
        # disk along with the input files
        with open(writer_metadata_path, "w") as outfile:
            json.dump(writer_metadata, outfile, indent=4, sort_keys=True)

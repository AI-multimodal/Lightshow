"""Contains all relevant classes and functions for constructing databases of
various materials and clusters. This includes pulling and processing data from
the `Materials Project Database <https://materialsproject.org/>`_, as well as
utilizing existing data the user may have on their hard drive."""

from datetime import datetime
import json
import os
from pathlib import Path
from shutil import copy2
from warnings import warn

from monty.json import MSONable
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester, MPRestError
from tqdm import tqdm

from lightshow import _get_API_key_from_environ
from lightshow import pymatgen_utils


def _delete_common_strings(old_list_of_strings):

    list_of_strings = [str(Path(xx).parent) for xx in old_list_of_strings]
    list_of_names = [str(Path(xx).name) for xx in old_list_of_strings]

    commonprefix = os.path.commonprefix(list_of_strings)
    new_p = [x[len(commonprefix) :] for x in list_of_strings]

    # Reverse every string in the list then do it again
    list_of_strings_reversed = [xx[::-1] for xx in new_p]
    commonsuffix = os.path.commonprefix(list_of_strings_reversed)
    new_p = [x[len(commonsuffix) :] for x in list_of_strings_reversed]

    final_p = [xx[::-1] for xx in new_p]

    return [
        str(Path(xx) / Path(name)) for xx, name in zip(final_p, list_of_names)
    ]


def _fetch_from_MP(mpr, mpid, metadata_keys):
    """Uses the provided MPID to fetch the structure data.

    Parameters
    ----------
    mpr : pymatgen.ext.matproj.MPRester
        Interface to the Materials Project REST API.
    mpid : str
        The specific mpid to pull.
    metadata_keys : list of str, optional
        The Materials Project metadata contains a huge amount of information.
        If not ``None``, these are the only keys that are kept in the pulled
        metadata object.

    Returns
    -------
    tuple
        The structure (:class:`pymatgen.core.structure.Structure`) of interest
        and the specified metadata.
    """

    metadata = mpr.get_doc(mpid)

    if metadata_keys is not None:
        metadata = {key: metadata[key] for key in metadata_keys}
    metadata["downloaded_at"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # The structure is precisely in the data pulled from get_doc:
    structure = Structure.from_dict(metadata.pop("structure"))

    return structure, metadata


def _from_mpids_list(mpids, api_key, metadata_keys, verbose=True):
    """Makes one large API call to the Materials Project database and pulls the
    relevant structural files given a list of Materials Project ID's (mpids).

    Parameters
    ----------
    mpids : list of str
        List of Materials Project IDs to pull.
    api_key : str, optional
        Materials Project API key.
    metadata_keys : bool, optional
        If True, will only
    verbose : bool, optional
        If True, will use tqdm to print a progress bar.

    Returns
    -------
    dict
        A dictionary containing the structures (with keys as the mpids) and
        metadata (with the same keys).
    """

    structures = dict()
    metadatas = dict()
    with MPRester(api_key) as mpr:
        for mpid in tqdm(mpids, disable=not verbose):
            try:
                structure, metadata = _fetch_from_MP(mpr, mpid, metadata_keys)
            except MPRestError as error:
                warn(f"MPRestError pulling mpid={mpid}, error: {error}")
                continue
            structures[mpid] = structure.get_primitive_structure()
            metadatas[mpid] = metadata

    return {"structures": structures, "metadata": metadatas}


def _get_api_key(api_key):
    if api_key is None:
        api_key = _get_API_key_from_environ()
    if api_key is None:
        raise ValueError(f"Invalid API key {api_key}")
    return api_key


class Database(MSONable):
    """Contains all materials and metadata for some database."""

    def cleanup_paths(self):
        """When loading data from disk, paths can become very repetitive. This
        method strips all common prefixes and suffixes from the keys of the
        structures and metadata properties."""

        old_keys = list(self._structures.keys())
        new_keys = _delete_common_strings(old_keys)
        for old_key, new_key in zip(old_keys, new_keys):
            self._structures[new_key] = self._structures.pop(old_key)
            self._metadata[new_key] = self._metadata.pop(old_key)

            if self._supercells_initialized:
                self._supercells[new_key] = self._supercells.pop(old_key)

    @classmethod
    def from_files(cls, root, filename="CONTCAR", cleanup_paths=True):
        """Searches for files matching the provided ``filename``, which can
        include wildcards, and assumes those files are structural files in CIF
        format. The names/ids of these files is given by the full directory
        structure where that file was found. For example, if
        ``root == "my_dir"``, ``filename == "CONTCAR"`` and we have a single
        structure file in ``my_dir/test/CONTCAR``, then the resulting
        structures will be ``{"my_dir/test/CONTCAR": struct}``. Similarly, if
        ``filename == "CONTCAR*", then all files of the form ``CONTCAR*`` will
        be found and used. The directory structure produced will be something
        like

        .. code-block:: python

            {
                "my_dir/test/CONTCAR1": struct1,
                "my_dir/test/CONTCAR2": struct2,
                ...
            }

        Parameters
        ----------
        root : str
            The directory in which to begin the search.
        filename : str, optional
            The files to search for. Uses ``rglob`` to recursively find any
            files matching ``filename`` within the provided directory.
        cleanup_paths : bool, optional
            If True, runs :class:`cleanup_paths()` after initializing the
            Database.

        Returns
        -------
        Database
        """

        structures = {
            str(Path(path.parent) / path.stem): Structure.from_file(
                path
            ).get_primitive_structure()
            for path in Path(root).rglob(filename)
        }

        # Check for any duplicate names
        names = list(structures.keys())
        names = [Path(name).name for name in names]
        if len(list(set(names))) < len(names):
            new_structures = dict()
            for key, value in structures.items():
                new_structures[str(Path(key).parent)] = value
            structures = new_structures

        metadata = {key: dict() for key in structures.keys()}

        kls = cls(structures=structures, metadata=metadata, supercells=dict())
        if cleanup_paths:
            kls.cleanup_paths()
        return kls

    @classmethod
    def from_materials_project(
        cls,
        query,
        query_type="mpids",
        api_key=None,
        metadata_keys=[
            "created_at",
            "blessed_tasks",
            "pseudo_potential",
            "spacegroup",
            "_id",
            "structure",
            "icsd_ids",
            "e_above_hull",
            "formation_energy_per_atom",
            "band_gap",
            "diel",
        ],
        verbose=True,
    ):
        """Constructs the :class:`.Database` object by pulling structures and
        metadata directly from the Materials Project. The following query types
        are allowed:

        * ``query_type == mpids``, the ``query`` argument is simply a list of
          Materials Project IDs.
        * ``query_type == patterns``, the ``query`` argument is a list
          of patterns, e.g. ``[Ti-O, Ti-O-*, Cu-O-*]``.
        * ``query_type == mp_query``, then the ``query`` argument is
          just a dict that is passed directly to the ``mpr.query`` method,
          where ``mpr`` is the MPRester object.

        Parameters
        ----------
        query : list or dict
            The query itself, the form of which is determined by
            ``query_type``.
        query_type : str, optional
            There are three different types of allowed queries: ``mpids``,
            ``patterns``, or ``mp_query``. See above for the allowed query
            types.
        api_key : str, optional
            Materials Project API key. If None (not provided), looks for the
            environment variable ``PMG_API_KEY``.
        metadata_keys : list, optional
            The Materials Project metadata contains a huge amount of
            information. If not ``None``, these are the only keys that are kept
            in the pulled metadata object.
        verbose : bool, optional
            If True, uses the :class:`tqdm` progress bar. Otherwise will
            silence it.

        Returns
        -------
        Database

        Examples
        --------
        Construct a :class:`.Database` via directly pulling certain materials
        by MPID.

        .. code-block:: python

            database = Database.from_materials_project(["mp-390", "mvc-1115"])

        Construct a :class:`.Database` via pulling all materials consistent
        with certain patterns. For example, to pull all binary and ternary
        titanium oxide compounds:

        .. code-block:: python

            database = Database.from_materials_project(
                ["Ti-O", "Ti-O-*"], query_type="patterns"
            )
        """

        api_key = _get_api_key(api_key)

        # Nothing to do here
        if query_type == "mpids":
            mpids = query

        # Convert the patterns into a list of mpids
        elif query_type == "patterns":
            mpids = []
            with MPRester(api_key) as mpr:
                for pattern in query:
                    data = mpr.get_data(pattern, data_type="vasp")
                    mpids.extend([xx["material_id"] for xx in data])

        # Convert the raw query itself to a list of mpids
        elif query_type == "mp_query":
            with MPRester(api_key) as mpr:
                materials_list = mpr.query(**query)
            mpids = [xx["material_id"] for xx in materials_list]

        # Otherwise we error and terminate
        else:
            raise ValueError(f"Unknown query {query_type}")

        # Get all of the data as a dictionary
        data = _from_mpids_list(
            mpids,
            api_key,
            metadata_keys,
            verbose=verbose,
        )

        return cls(
            structures=data["structures"],
            metadata=data["metadata"],
            supercells=dict(),
        )

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
        metadata=dict(),
        supercells=dict(),
        supercell_cutoff=None,
        inequivalent_sites_initialized=False,
        supercells_initialized=False,
    ):
        """Initializer for the Database class. Note it is recommended to use
        the classmethods to initialize this object."""

        self._structures = structures
        self._metadata = metadata
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

    def write(
        self,
        root,
        absorbing_atoms=None,
        options=[],
        pbar=True,
        copy_script=None,
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
            performed.
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
        """

        self._setup_preliminary_attributes()

        root = str(Path(root).resolve())  # Get absolute path

        if not isinstance(absorbing_atoms, list):
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

            for absorbing_atom in absorbing_atoms:

                primitive_info = self._metadata[key]["primitive"]
                supercell_info = self._metadata[key]["supercell"]

                # If inequiv is None, that means that the absorbing_atom was not
                # specified (absorbing_atom is None)
                inequiv = self._get_site_indexes_matching_atom(
                    primitive_info, absorbing_atom
                )
                inequiv_sc = self._get_site_indexes_matching_atom(
                    supercell_info, absorbing_atom
                )

                if inequiv_sc is not None and inequiv is not None:
                    index_mapping = {k: v for k, v in zip(inequiv, inequiv_sc)}
                else:
                    index_mapping = None

                # If the for VASP and XSpectra calculations, use supercell;
                # otherwise, use unit cell structure
                # test if band_gap and diel in the self._metadata[key].keys()
                # if yes, read the bandgap and diel for OCEAN
                # if no, ignore them
                kwargs = {
                    "structure_sc": supercell,
                    "structure_uc": self._structures[key],
                    "sites": inequiv,
                    "index_mapping": index_mapping,
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

        # Save a metadata file (not a serialized version of this class) to
        # disk along with the input files
        with open(writer_metadata_path, "w") as outfile:
            json.dump(writer_metadata, outfile, indent=4, sort_keys=True)

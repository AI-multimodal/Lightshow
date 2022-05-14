"""Contains all relevant classes and functions for constructing databases of
various materials and clusters. This includes pulling and processing data from
the `Materials Project Database <https://materialsproject.org/>`_, as well as
utilizing existing data the user may have on their hard drive."""

from datetime import datetime
from pathlib import Path

from monty.json import MSONable
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester, MPRestError
from tqdm import tqdm

from lightshow import _get_API_key_from_environ


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


def _from_mpids_list(
    mpids, api_key, metadata_keys, verbose=True, suppress_MPRestError=False
):
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
    suppress_MPRestError : bool, optional
        If True, will ignore any ``MPRestError`` thrown during the process of
        pulling data.

    Returns
    -------
    dict
        A dictionary containing the structures (with keys as the mpids)
        metadata (with the same keys), and any mpids that were missed due to a
        :class:`pymatgen.ext.matproj.MPRestError`.
    """

    structures = dict()
    metadatas = dict()
    errors = []
    with MPRester(api_key) as mpr:
        for mpid in tqdm(mpids, disable=not verbose):
            try:
                structure, metadata = _fetch_from_MP(mpr, mpid, metadata_keys)
            except MPRestError:
                errors.append(mpid)
                continue
            structures[mpid] = structure
            metadatas[mpid] = metadata

    return {"structures": structures, "metadata": metadatas, "errors": errors}


def _get_api_key(api_key):
    if api_key is None:
        api_key = _get_API_key_from_environ()
    if api_key is None:
        raise ValueError(f"Invalid API key {api_key}")
    return api_key


class Database(MSONable):
    """Contains all materials and metadata for some database."""

    @classmethod
    def from_files(cls, root, filename="CONTCAR", api_key=None, verbose=True):
        """Searches for files matching the provided ``filename``, and assumes
        those files are structural files in CIF format. The names/ids of these
        files is given by the full directory structure where that file was
        found. For example, if ``root == "my_dir"``, ``filename == "CONTCAR"``
        and we have a single structure file in ``my_dir/test/CONTCAR``, then
        the resulting structures will be ``{"my_dir/test": struct}``.

        Parameters
        ----------
        root : str
            The directory in which to begin the search.
        filename : str, optional
            The files to search for. Uses ``rglob`` to recursively find any
            files matching ``filename`` within the provided directory.
        api_key : str, optional
            Materials Project API key. If None (not provided), looks for the
            environment variable ``PMG_API_KEY``.
        verbose : bool, optional
            If True, will use tqdm to print a progress bar.

        Returns
        -------
        Database
        """

        api_key = _get_api_key(api_key)
        structures = {
            str(path.parent): Structure.from_file(path)
            for path in Path(root).rglob(filename)
        }
        return cls(structures, None)

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
        data = _from_mpids_list(mpids, api_key, metadata_keys, verbose=verbose)

        # Instantiate the class and return
        return cls(data["structures"], data["metadata"])

    @property
    def structures(self):
        """A dictionary of :class:`pymatgen.core.structure.Structure` objects.
        The keys are the IDs of the structure. In the case of data pulled from
        the Materials Project, these are the MPIDs. Otherwise, they are simply
        strings encoding some information about the origins of the structures.

        Returns
        -------
        dict
        """

        return self._structures

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

    def __init__(self, structures, metadata):
        """Initializer for the Database class. Note it is recommended to use
        the classmethods to initialize this object."""

        self._structures = structures
        self._metadata = metadata

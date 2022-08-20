==========
Quickstart
==========

.. note::

    In addition to the quickstart guide below, we also have tutorial notebooks under the ``notebooks`` directory in the Lightshow GitHub repository.

The database
============

Everything in **Lightshow** is built around materials databases, specifically the ``lightshow.database.Database`` class. This object is a lightweight container for Pymatgen ``Structure`` objects, and can be instantiated in a few different ways.

Begin by importing the ``Database`` object.

.. code-block:: python

    from lightshow import Database

.. tip::

    Almost every object in **Lightshow** can be json-serialized using the ``MSONable`` syntax. For example,

        .. code-block:: python

            d = database.as_dict()
            new_database = Database.from_dict(d)
            my_json_rep = database.to_json()

            import json
            with open("my_database.json", "w") as f:
                json.dump(my_json_rep, f)

Databases from the Materials Project
------------------------------------

The database is designed to be constructed via classmethods. The primary classmethod we recommend using is ``from_materials_project``. It interfaces directly with the ``MPRester`` object to make queries and pull data locally. There are three ways to do this:

#. Directly pull a list of Materials Project IDs
   
   .. code-block:: python

        database = Database.from_materials_project(
            query=["mp-390", "mvc-11115"],
            query_type="mpids",
            api_key=API_KEY
        )

#. Use the wildcard/patterns syntax to pull larger databases matching certain patterns. For example, to pull all binary and ternary Ti-O structures,
   
   .. code-block:: python

       database = Database.from_materials_project(
            query=["Ti-O", "Ti-O-*"],
            query_type="patterns",
            api_key=API_KEY
        )

#. Via direct REST query. See the appropriate Pymatgen docs `here <https://pymatgen.org/pymatgen.ext.matproj.html?highlight=mprester#pymatgen.ext.matproj.MPRester.query>`__ for more details.

   .. code-block:: python

       database = Database.from_materials_project(
            query={
                "criteria": ...,
                "properties": ...,
                ...
            },
            query_type="mp_query",
            api_key=API_KEY
        )

.. note::

    While the `Pymatgen API Key <https://legacy.materialsproject.org/open>`_ can be provided manually during the use of ``from_materials_project``, we highly recommend setting it in the environment variable ``PMG_API_KEY``. If ``api_key`` in the above arguments is not provided or is set to ``None``, **Lightshow** will check for this environment variable and use it instead.

Once the ``database`` has been built, three properties are accessible:

#. ``database.structures`` is a dictionary of Materials Project IDs (MPIDS) as keys and the ``Structure`` as values.
#. ``database.metadata`` is a dictionary of MPIDs as keys and dictionaries containing metadata as values.
#. ``database.errors`` is a dictionary containing any errors logged during the usage of the ``Database`` object. 

.. note::

    We fully document all "public" functions, classes and methods. Documentation can be easily accessed through the Lightshow API reference (see the sidebar) or by doing e.g. ``Database.from_materials_project?`` in a Jupyter Notebook or iPython instance.


Databases from disk
-------------------
It is also possible to construct the ``Database`` from data on disk. This method will not fill the ``metadata`` property, though, which might force **Lightshow** to rely on default parameter values for certain types of input files.

    .. code-block:: python

        database = Database.from_files(
            root="search/here/for/files",
            filename="CONTCAR"
        )

The code above will look recursively in the provided directory for files matching the ``filename`` argument, and will attempt to read those into a ``Structure`` object via Pymatgen's ``Structure.from_file`` classmethod. The keys to the ``database.structures`` property will be the path to the parent directory containing the structure file instead of the MPID.


Parameters
==========

Our primary common abstraction is that of the spectroscopy simulation parameters. These control every aspect of the input files to be written and are specific to each type of code. However, while all options are exposed for the user, sensible defaults are also provided, making it straightforward to get started. Currently, we provide support for 5 different codes: FEFF, VASP, EXCITING, OCEAN and Xspectra.

You can begin by importing the simulation code-specific parameter objects,

.. code-block:: python

    from lightshow import (
        FEFFParameters,
        VASPParameters,
        OCEANParameters,
        XSpectraParameters,
        EXCITINGParameters
    )

which we will go over one-by-one.

.. _feff-reference:

FEFF
----

.. note::

    See `here <https://feff.phys.washington.edu/feff/Docs/feff9/feff90/feff90_users_guide.pdf>`__ for the FEFF9 documentation.

There are three primary arguments for the ``FEFFParameters`` object: the ``cards``, ``edge`` and ``radius``. For example,

.. code-block:: python

    feff_params = FEFFParameters(
        cards={
            "S02": "0",
            "COREHOLE": "RPA",
            "CONTROL": "1 1 1 1 1 1",
            "XANES": "4 0.04 0.1",
            "SCF": "7.0 0 100 0.2 3",
            "FMS": "9.0 0",
            "EXCHANGE": "0 0.0 0.0 2",
            "RPATH": "-1"
        },
        edge="K",
        radius=10.0
    )

``cards`` is a catch-all input which is written directly to the preamble of the ``feff.inp`` file. Essentially, any parameter can be provided here, and should be provided as strings (both keys and values). A complete list of allowed "control cards" can be found on page 69 of the `FEFF9 documentation <https://feff.phys.washington.edu/feff/Docs/feff9/feff90/feff90_users_guide.pdf>`__. Note that certain cards, while required, are not directly passed using ``cards`` above. For example, the ``POTENTIALS`` card is automatically written.

``edge`` determines the x-ray absorption edge of the calculation. Particulars are noted on page 89 of the `FEFF9 documentation <https://feff.phys.washington.edu/feff/Docs/feff9/feff90/feff90_users_guide.pdf>`__.

.. warning::

    According to the FEFF9 documentation, M-shells or higher are not well tested. **Lightshow** will provide a warning if the user sets these edges.

``radius`` is a critical parameter that sets the cluster size. For each absorbing atom, a radius of ``radius`` Å is taken around that absorbing atom, a supercell is appropriately constructed, and then truncated such that the only atoms contained in the ``feff.inp`` file are at most ``radius`` Å away from the absorber. Note that in this sense, ``radius`` controls much of the computational expense of the FEFF calculation.

The remainder of the ``feff.inp`` file is constructed automatically, and to some degree leverages `Pymatgen's FEFF IO module <https://pymatgen.org/pymatgen.io.feff.inputs.html>`__.


VASP
----

The added complexity of the VASP input files necessitates slightly more complicated syntax on the side of **Lightshow**. During any VASP run, there are four objects that are required in the working directory before running the VASP executable: INCAR, KPONTS, POTCAR and POSCAR, representing the general input file parameters, k-points parameters, pseudopotential files, and structure files, respectively.

.. note::

    The VASP documentation can be found `here <https://www.vasp.at/wiki/index.php/The_VASP_Manual>`__.

The general ``VASPParameter`` object structure looks something like this:

.. code-block:: python

    vasp_parameters = VASPParameters(
        incar=...,
        edge="K",
        potcar_directory=None
    )

where for now we have suppressed some sensible defaults which are discussed later. The primary information required to instantiate the :class:`lightshow.parameters.vasp.VASPParameters` object are the ``incar``, ``edge``, and ``potcar_directory`` arguments.

INCAR sets the parameters for the INCAR input file. It can either take a Python dictionary, or :class:`lightshow.parameters.vasp.Incar` object. The only parameter that can be overwritten is ``incar["NBANDS"]``. If this INCAR parameter is ``None``, **Lightshow** will attempt to use a default method to estimate a good number of bands for the VASP calculation. This is discussed more in :ref:`customize-n-bands` below. We provide sensible defaults for the INCAR files in :class:`lightshow.parameters.vasp`.

``edge`` sets the x-ray absorption edge. See :ref:`feff-reference`.


Pseudopotentials
^^^^^^^^^^^^^^^^

``potcar_directory`` points **Lightshow** to a directory containing VASP pseudopotential files. The handling of these files can be confusing, hence we outline how **Lightshow** handles them here in detail.

.. warning::

    VASP POTCAR (potential) files are under the VASP license and thus are not included in **Lightshow**. In order to use VASP and the potential files, you must have a VASP license. See the `VASP Website <https://www.vasp.at>`__ for more details.

The :class:`lightshow.parameters.vasp.PotcarConstructor` handles creating the POTCAR file when writing the input files. The default parameters of this object can be overwritten through various other arguments in :class:`lightshow.parameters.vasp.VASPParameters`, but the defaults are recommended.

.. tip::

    If ``potcar_directory`` is ``None``, **Lightshow** will attempt to read this from an environment variable ``VASP_POTCAR_DIRECTORY``.

The directory that ``potcar_directory`` points to should contain files of the form of the values in :class:`lightshow.parameters.vasp.VASP_POTCAR_DEFAULT_ELEMENT_MAPPING`. Specific values for these mappings, which map element types to specific potential files in the directory provided, can be overwritten by setting ``potcar_element_mapping``. These provided values will only override the keys provided.



OCEAN
-----

TODO


EXCITING
--------
.. note::

    See `here <http://exciting.wikidot.com/ref:input>`__ for the EXCITING documentation.

There are three required primary arguments for the ``EXCITINGParameters`` object: the ``cards``, ``species_directory`` and ``edge``. For example, the general ``VASPParameter`` object structure looks something like this:

.. code-block:: python

    exciting_params = EXCITINGParameters(
        cards={
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
                    "gqmax": 4.0,
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
             
        },
        species_directory = 'this/is/not/a/directory'
        edge="K",
    )


``cards`` is a catch-all input which is written directly to the preamble of the ``input.xml`` file. Essentially, any parameter can be provided here, and should be provided as strings (both keys and values). Note that certain cards, while required, are not directly passed using ``cards`` above. An example is the ``species_directory`` discussed below.

``species_directory`` lets the user to specify where the species files are located. Usually, one can find these files in the exciting souce code under the ``species`` directory. If ``species_directory`` is not set, a warning will show up indicating the users should copy the corresponding species files to the working directory, e.g. where the ``input.xml`` file is generated.

``edge`` sets the x-ray absorption edge. See :ref:`feff-reference`.


XSpectra
--------

TODO


Advanced
========

.. _customize-k-points:

Customize k-points
------------------

.. _customize-n-bands:

Customize the number of bands
-----------------------------





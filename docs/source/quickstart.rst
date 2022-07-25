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

#. Via direct REST query. See the appropriate Pymatgen docs [#pmgREST]_ for more details.

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

FEFF
----

.. note::

    See here [#feff9docs]_ for the FEFF9 documentation.

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

``cards`` is a catch-all input which is written directly to the preamble of the ``feff.inp`` file. Essentially, any parameter can be provided here, and should be provided as strings (both keys and values). A complete list of allowed "control cards" can be found on page 69 of the FEFF9 documentation [#feff9docs]_. Note that certain cards, while required, are not directly passed using ``cards`` above. For example, the ``POTENTIALS`` card is automatically written.

``edge`` determines the x-ray absorption edge of the calculation. Particulars are noted on page 89 of the FEFF9 documentation [#feff9docs]_. 

.. warning::

    According to the FEFF9 documentation, M-shells or higher are not well tested. **Lightshow** will provide a warning if the user sets these edges.

``radius`` is a critical parameter that sets the cluster size. For each absorbing atom, a radius of ``radius`` Å is taken around that absorbing atom, a supercell is appropriately constructed, and then truncated such that the only atoms contained in the ``feff.inp`` file are at most ``radius`` Å away from the absorber. Note that in this sense, ``radius`` controls much of the computational expense of the FEFF calculation.

The remainder of the ``feff.inp`` file is constructed automatically, and to some degree leverages Pymatgen's ``feff.io`` module [#pmgFEFFio]_.


VASP
----

TODO


OCEAN
-----

TODO


EXCITING
--------

TODO


XSpectra
--------

TODO

----

.. [#pmgREST] Pymatgen's REST query docs can be found at `pymatgen.org/pymatgen.ext.matproj.html?highlight=mprester#pymatgen.ext.matproj.MPRester.query <https://pymatgen.org/pymatgen.ext.matproj.html?highlight=mprester#pymatgen.ext.matproj.MPRester.query>`_
.. [#feff9docs] FEFF9 documentation can be found at `feff.phys.washington.edu/feff/Docs/feff9/feff90/feff90_users_guide.pdf <https://feff.phys.washington.edu/feff/Docs/feff9/feff90/feff90_users_guide.pdf>`_.
.. [#pmgFEFFio] Pymatgen's ``feff.io.inputs`` module can be found at `pymatgen.org/pymatgen.io.feff.inputs.html <https://pymatgen.org/pymatgen.io.feff.inputs.html>`_.




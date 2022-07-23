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

#. Via direct REST query. See `here <https://pymatgen.org/pymatgen.ext.matproj.html?highlight=mprester#pymatgen.ext.matproj.MPRester.query>`_ for more details.

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

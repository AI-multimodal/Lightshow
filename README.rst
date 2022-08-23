Lightshow
=========

.. inclusion-marker-LIGHTSHOW-begin

**Lightshow** is a Python library meant for easily generating computational spectroscopy input files.

Often, it can be a daunting task to create comprehensive, well documented databases of materials structures and their x-ray absorption spectra. **Lightshow** solves this problem, allowing new users to choose sensible defaults for their calculations, while simultaneously exposing all functionality for experts.

**Lightshow** aims to provide a "one-stop-shop" for input file generation, and currently supports the following codes:

- FEFF
- VASP
- OCEAN
- EXCITING
- Xspectra

with more on the way! The software is intended to be user-friendly, completely documented and tested, and extendable for those users who wish to add additional spectroscopy functionalities. There are also a few comprehensive tutorials to help you get started.

.. inclusion-marker-LIGHTSHOW-end

Installation
------------

.. inclusion-marker-LIGHTSHOW-installation-begin

Users
^^^^^
To simply use the software, install it as you would any Python package: `pip install lightshow`. **COMING SOON!**

Developers
^^^^^^^^^^
If you wish to help us improve **Lightshow**, you should fork a copy of our repository, clone to your local machine, and then proceed with setting up the following:

Create and activate a fresh virtual environment, e.g.

.. code-block:: bash
    
    conda create -n py3.9 python=3.9 && conda activate py3.9

It is highly recommended that you also install the pre-commit hooks. This will help you avoid failing the black and flake8 tests that are required as part of our CI testing suite.

.. code-block:: bash

    pre-commit install

`Install the development requirements <https://github.com/pypa/pip/issues/8049#issuecomment-633845028>`_. We use the helper script ``build.sh`` to parse the ``pyproject.toml`` file and only install the specified packages needed for development (note this does not actually install **Lightshow**).

.. code-block:: bash
    
    bash build.sh install-dev-requirements



.. - Setup the pre-commit hooks ``pre-commit install``


.. inclusion-marker-LIGHTSHOW-installation-end

.. inclusion-marker-LIGHTSHOW-funding-begin

Funding acknowledgement
-----------------------
This research is based upon work supported by the U.S. Department of Energy, Office of Science, Office Basic Energy Sciences, under Award Number FWP PS-030. This research used resources of the Center for Functional Nanomaterials (CFN), which is a U.S. Department of Energy Office of Science User Facility, at Brookhaven National Laboratory under Contract No. DE-SC0012704.

.. inclusion-marker-LIGHTSHOW-funding-end

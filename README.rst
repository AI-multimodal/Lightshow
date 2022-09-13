.. inclusion-marker-LIGHTSHOW-begin

.. image:: https://raw.githubusercontent.com/AI-multimodal/Lightshow/master/docs/_static/images/lightshow.jpg
    :align: right
    :alt: sysfs line plot
    :width: 1000px

.. image:: https://github.com/AI-multimodal/Lightshow/actions/workflows/ci.yml/badge.svg
 :target: https://github.com/AI-multimodal/Lightshow/actions/workflows/ci.yml

.. image:: https://codecov.io/gh/AI-multimodal/Lightshow/branch/master/graph/badge.svg?token=CW7BMFA5O7 
 :target: https://codecov.io/gh/AI-multimodal/Lightshow
 
.. image:: https://app.codacy.com/project/badge/Grade/d31a4e18672c4d71bbaafa719181c140
 :target: https://www.codacy.com/gh/AI-multimodal/Lightshow/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=AI-multimodal/Lightshow&amp;utm_campaign=Badge_Grade

=============

**Lightshow** is a Python library for easily generating computational spectroscopy input files.

.. inclusion-marker-LIGHTSHOW-end

.. inclusion-marker-LIGHTSHOW-overview-begin

Often, it can be a daunting task to create comprehensive, well documented databases of materials structures and their x-ray absorption spectra. **Lightshow** solves this problem, allowing new users to choose sensible defaults for their calculations, while simultaneously exposing all functionality for experts.

**Lightshow** aims to provide a "one-stop-shop" for input file generation, and currently supports the following codes:

- FEFF
- VASP
- OCEAN
- EXCITING
- Xspectra

with more on the way! The software is intended to be user-friendly, completely documented and tested, and extendable for those users who wish to add additional spectroscopy functionalities. There are also a few comprehensive tutorials to help you get started.

.. inclusion-marker-LIGHTSHOW-overview-end


Installation
------------

.. inclusion-marker-LIGHTSHOW-installation-begin

Users
^^^^^
To simply use the software, install it as you would any Python package: 

.. code-block:: bash

    pip install lightshow


Developers
^^^^^^^^^^
If you wish to help us improve **Lightshow**, you should fork a copy of our repository, clone to your local machine, and then proceed with setting up the following:

Create and activate a fresh virtual environment, e.g.

.. code-block:: bash
    
    conda create -n py3.9 python=3.9 && conda activate py3.9

It is highly recommended that you also install the pre-commit hooks. This will help you avoid failing the black and flake8 tests that are required as part of our CI testing suite.

.. code-block:: bash

    pre-commit install

We use helper scripts to parse the ``pyproject.toml`` file and install only specific packages required for certain parts of development. For development, we recommend installing all dependencies:

.. code-block:: bash
    
    bash scripts/install.sh       # Install Lightshow's core dependencies
    bash scripts/install.sh test  # Install the test requirements only
    bash scripts/install.sh doc   # Install requirements for building the docs

.. inclusion-marker-LIGHTSHOW-installation-end

.. inclusion-marker-LIGHTSHOW-funding-begin

Funding acknowledgement
-----------------------
This research is based upon work supported by the U.S. Department of Energy, Office of Science, Office Basic Energy Sciences, under Award Number FWP PS-030. This research used resources of the Center for Functional Nanomaterials (CFN), which is a U.S. Department of Energy Office of Science User Facility, at Brookhaven National Laboratory under Contract No. DE-SC0012704.

.. inclusion-marker-LIGHTSHOW-funding-end

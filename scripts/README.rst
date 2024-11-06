Helper scripts
==============

This directory contains a variety of scripts which are used not only by the continuous integration (CI) pipeline, but can also be used by contributors during the development process.


Installation
------------
`install.sh` is a helper script for parsing the `pyproject.toml` file (see `here <https://github.com/pypa/pip/issues/8049>`_ for the inspiration for these scripts) and installing the required packages as specified there. As also explained in the `README <../README.rst>`_, we can install all or part of Lightshow's dependencies via the following scripts:

.. code-block:: bash
    
    bash scripts/install.sh       # Install Lightshow's core dependencies
    bash scripts/install.sh test  # Install the test requirements only
    bash scripts/install.sh doc   # Install requirements for building the docs

Building
--------
We provide two self-contained scripts for building the Lightshow source and docs. These scripts are `build_project.sh` and `build_docs.sh`, respectively.


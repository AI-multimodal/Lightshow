============
Installation
============

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

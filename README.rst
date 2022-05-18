Lightshow
=========

.. inclusion-marker-LIGHTSHOW-begin

``Lightshow`` is a one-stop-shop for generating computational spectroscopy input files. ``Lightshow`` heavily wraps Pymatgen, is based largely on the familiar ``Structure`` object, and adheres to the following design philosophies:

- **User friendliness:** All input parameters for all spectroscopies are exposed through the API, but the key use cases are highlighted and brought to the forefront.
- **Complete documentation and tutorials:** Essentially every use case will come with its own set of tutorials and detailed instructions. It goes without saying, all code is documented.
- **Extendability at the forefront:** We don't plan on stopping at the spectroscopy methods that already exist in the code, and the software is modularized such that adding new computational methods should be straightforward.

.. inclusion-marker-LIGHTSHOW-end

Installation
------------

Users
^^^^^
To simply use the software, install it as you would any Python package: `pip install lightshow`. **COMING SOON!**

Developers
^^^^^^^^^^
If you wish to help us improve Lightshow, you should fork a copy of our repository, clone to disk, and then proceed with setting up the following:

- Create a fresh virtual environment, e.g. ``conda create -n py3.9 python=3.9``.
- Install the development requirements, ``pip install -r requirements-dev.txt``
- Setup the pre-commit hooks ``pre-commit install``
- If you want to install the package to your default paths, you can do this in "developer mode" by running ``pip install -e ".[dev]"``


.. inclusion-marker-LIGHTSHOW-funding-begin

Funding acknowledgement
-----------------------
This research is based upon work supported by the U.S. Department of Energy, Office of Science, Office Basic Energy Sciences, under Award Number FWP PS-030.

.. inclusion-marker-LIGHTSHOW-funding-end

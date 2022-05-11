# Lightshow (temporary title?)
A one-stop-shop code for writing input files for computational spectroscopy calculations. Our software has the following key features:
* Consistent, seamless API
* Everything is serializable, allowing for maximal reproducibility
* Everything is logged for the user's convenience
* We interface primarily with the Materials Project, and utilize their native `Structure` object whenever possible, allowing for maximal transferability and extendability

# Install

## Users
To simply use the software, install it as you would any Python package: `pip install lightshow`. **COMING SOON!**

## Developers
If you wish to help us improve Lightshow, you should fork a copy of our repository, clone to disk, and then proceed with setting up the following:

* Create a fresh virtual environment, e.g. `conda create -n py3.9 python=3.9`.
* Install the development requirements, `pip install -r requirements-dev.txt`
* Setup the pre-commit hooks `pre-commit install`
* If you want to install the package to your default paths, you can do this in "developer mode" by running `pip install -e ".[dev]"`

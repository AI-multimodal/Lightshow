<div align="center">
    
<!-- ![sysfs line plot](https://raw.githubusercontent.com/AI-multimodal/Lightshow/master/docs/_static/images/lightshow.jpg) -->

# Lightshow

<!--badges-start-->
[![image](https://joss.theoj.org/papers/a9cabcd7f4b85a926a797997c6622b43/status.svg)](https://joss.theoj.org/papers/a9cabcd7f4b85a926a797997c6622b43)
[![image](https://github.com/AI-multimodal/Lightshow/actions/workflows/ci.yml/badge.svg)](https://github.com/AI-multimodal/Lightshow/actions/workflows/ci.yml)
[![image](https://codecov.io/gh/AI-multimodal/Lightshow/branch/master/graph/badge.svg?token=CW7BMFA5O7)](https://codecov.io/gh/AI-multimodal/Lightshow)
[![image](https://app.codacy.com/project/badge/Grade/d31a4e18672c4d71bbaafa719181c140)](https://www.codacy.com/gh/AI-multimodal/Lightshow/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=AI-multimodal/Lightshow&amp;utm_campaign=Badge_Grade) <br>
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![python](https://img.shields.io/badge/-Python_3.9+-blue?logo=python&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Downloads](https://static.pepy.tech/badge/lightshow)](https://pepy.tech/project/lightshow)
<!--badges-end-->
</div>

    
------------------------------------------------------------------------

<!--lightshow-intro-start-->
**Lightshow** is a Python library for easily generating computational
spectroscopy input files. If you use our code, please consider citing our [manuscript](https://doi.org/10.21105/joss.05182) in the Journal of Open Source Software and our multi-code benchmark [paper](https://doi.org/10.1103/PhysRevMaterials.8.013801), which provides the methods and default parameters used in Lightshow.

Often, it can be daunting to create comprehensive, well documented
databases of materials structures and their x-ray absorption spectra.
**Lightshow** solves this problem, allowing new users to choose sensible
defaults for their calculations, while simultaneously exposing all
functionality for experts.

**Lightshow** aims to provide a \"one-stop-shop\" for input file
generation, and currently supports the following codes:

-   [FEFF](https://feff.phys.washington.edu)
-   [VASP](https://www.vasp.at)
-   [OCEAN](https://www.nist.gov/services-resources/software/ocean)
-   [EXCITING](https://exciting-code.org)
-   [Xspectra](https://gitlab.com/QEF/q-e/-/tree/master/XSpectra)

with more on the way! The software is intended to be user-friendly,
extensively documented and tested, and extendable for those users who
wish to add additional spectroscopy functionalities. There are also a
few comprehensive tutorials to help you get started.

<!--lightshow-intro-end-->


## Tutorials

We offer a few tutorials to get you started (with more on the way!)

| # | Tutorial | Description | Notebook |
| ------------- | ------------- | ------------- | ------------- |
| 1 | Basic usage | Pull structure data from Materials Project, write input files for all codes | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AI-multimodal/Lightshow/blob/master/notebooks/00_basic_usage.ipynb) |
| 2 | Spectral broadening | Utilities for broadening spectra |[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AI-multimodal/Lightshow/blob/master/notebooks/01_Ti_K_anatase_broaden.ipynb) |


## Installation

<!--standard-installation-start-->

To install Lightshow, simply use `pip`

``` bash
pip install lightshow
```

Make sure you've set your Materials Project v2 API key as well! You can find some documentation on how to query data [here](https://docs.materialsproject.org/downloading-data/using-the-api/querying-data) and how to set up your new API key [here](https://next-gen.materialsproject.org/api).

```bash
export MP_API_KEY="your_Materials_Project_v2_API_key"
```

(or preferably, add `MP_API_KEY` to your bash profile).

More details can be found at our [documentation](https://ai-multimodal.github.io/Lightshow/installation.html).

<!--standard-installation-end-->


### Development installation


<!--dev-installation-start-->

For developers: after cloning Lightshow locally, install `pre-commit` via

```bash
pip install pre-commit
pre-commit
pre-commit install
```


and check that the tests (below) work correctly (these can be run via `pytest`). After cloning, simply use
```bash
pytest lightshow/_tests
```

or with coverage

```bash
pytest -v --cov --cov-report xml lightshow/_tests
```

We use helper scripts to parse the ``pyproject.toml`` file and install only specific packages required for certain parts of development. For development, we recommend installing all dependencies:

```bash
bash scripts/install.sh       # Install Lightshow's core dependencies
bash scripts/install.sh test  # Install the test requirements only
bash scripts/install.sh doc   # Install requirements for building the docs
```


<!--dev-installation-end-->


## Contributing 

We welcome any and all contributions by the community, including pull requests, bug reports, etc.
Please see our [contributing](https://github.com/AI-multimodal/Lightshow/blob/master/CONTRIBUTING.md) document for more details!

### Adding new spectroscopy codes

Adding new spectroscopy codes requires one to inherit the `_BaseParameters` class from `lightshow.parameters._base`.
The new `Parameters(_BaseParameters)` object should have a `write()` method, which must take a target directory as an
argument, as well as any other keyword arguments require to write the input file (most notably, the Pymatgen structure,
either a unit or super cell). The `name` property must also be defined (corresponding to the name of the calculation,
e.g. "VASP"). Take a look at the examples in the `lightshow.parameters` module for more details!


## Funding acknowledgement

<!--funding-start-->
This research is based upon work supported by the U.S. Department of
Energy, Office of Science, Office Basic Energy Sciences, under Award
Number FWP PS-030. This research used resources of the Center for
Functional Nanomaterials (CFN), which is a U.S. Department of Energy
Office of Science User Facility, at Brookhaven National Laboratory under
Contract No. DE-SC0012704. This work received partial funding by the German Research Foundation
(DFG) through the CRC 1404 (FONDA), Projektnummer 414984028, and the NFDI consortium FAIRmat – project 460197019.
C.V. acknowledges support by the Department of Energy, Basic Energy Sciences, Materials Science and Engineering Division,
through the Midwest Integrated Center for Computational Materials (MICCoM).

The Software resulted from work developed under a U.S. Government
Contract No. DE-SC0012704 and are subject to the following terms: the
U.S. Government is granted for itself and others acting on its behalf a
paid-up, nonexclusive, irrevocable worldwide license in this computer
software and data to reproduce, prepare derivative works, and perform
publicly and display publicly.

THE SOFTWARE IS SUPPLIED \"AS IS\" WITHOUT WARRANTY OF ANY KIND. THE
UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND THEIR
EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, TITLE OR NON-INFRINGEMENT, (2) DO NOT ASSUME
ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF THE
SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4) DO NOT WARRANT
THAT THE SOFTWARE WILL FUNCTION UNINTERRUPTED, THAT IT IS ERROR-FREE OR
THAT ANY ERRORS WILL BE CORRECTED.

IN NO EVENT SHALL THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
ENERGY, OR THEIR EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF ANY KIND OR
NATURE RESULTING FROM EXERCISE OF THIS LICENSE AGREEMENT OR THE USE OF
THE SOFTWARE.
<!--funding-end-->

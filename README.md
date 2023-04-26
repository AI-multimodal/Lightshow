<div align="center">

![sysfs line plot](https://raw.githubusercontent.com/AI-multimodal/Lightshow/master/docs/_static/images/lightshow.jpg)

[![image](https://github.com/AI-multimodal/Lightshow/actions/workflows/ci.yml/badge.svg)](https://github.com/AI-multimodal/Lightshow/actions/workflows/ci.yml)
[![image](https://codecov.io/gh/AI-multimodal/Lightshow/branch/master/graph/badge.svg?token=CW7BMFA5O7)](https://codecov.io/gh/AI-multimodal/Lightshow)
[![image](https://app.codacy.com/project/badge/Grade/d31a4e18672c4d71bbaafa719181c140)](https://www.codacy.com/gh/AI-multimodal/Lightshow/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=AI-multimodal/Lightshow&amp;utm_campaign=Badge_Grade)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![python](https://img.shields.io/badge/-Python_>=3.7-blue?logo=python&logoColor=white)](https://github.com/pre-commit/pre-commit) <br>
[![image](https://joss.theoj.org/papers/a9cabcd7f4b85a926a797997c6622b43/status.svg)](https://joss.theoj.org/papers/a9cabcd7f4b85a926a797997c6622b43)
[![image](https://zenodo.org/badge/DOI/10.48550/arXiv.2211.04452.svg)](https://doi.org/10.48550/arXiv.2211.04452)

</div>
    
------------------------------------------------------------------------

**Lightshow** is a Python library for easily generating computational
spectroscopy input files.

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

# Tutorials

We offer a few tutorials to get you started (with more on the way!)

| # | Tutorial | Description | Notebook |
| ------------- | ------------- | ------------- | ------------- |
| 1 | Basic usage | Pull structure data from Materials Project, write input files for all codes | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AI-multimodal/Lightshow/blob/master/notebooks/00_basic_usage.ipynb) |
| 2 | Spectral broadening | Utilities for broadening spectra |[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/AI-multimodal/Lightshow/blob/master/notebooks/01_Ti_K_anatase_broaden.ipynb) |


# Installation

To install Lightshow, simply use `pip`

``` bash
pip install lightshow
```

More details can be found at our [documentation](https://ai-multimodal.github.io/Lightshow/installation.html).


# Funding acknowledgement

This research is based upon work supported by the U.S. Department of
Energy, Office of Science, Office Basic Energy Sciences, under Award
Number FWP PS-030. This research used resources of the Center for
Functional Nanomaterials (CFN), which is a U.S. Department of Energy
Office of Science User Facility, at Brookhaven National Laboratory under
Contract No. DE-SC0012704. This work received partial funding by the German Research Foundation
(DFG) through the CRC 1404 (FONDA), Projektnummer 414984028, and the NFDI consortium FAIRmat – project 460197019.
C.V. acknowledges support by the Department of Energy, Basic Energy Sciences, Materials Science and Engineering Division,
through the Midwest Integrated Center for Computational Materials (MICCoM).

## Disclaimer

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

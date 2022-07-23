---
title: 'Lightshow: a Python package for writing computational x-ray absorption spectroscopy input files'
tags:
  - Python
  - computational spectroscopy
  - condensed matter physics
  - materials science
authors:
  - name: Matthew R. Carbone
    orcid: 0000-0002-5181-9513
    equal-contrib: true
    corresponding: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Fanchen Meng
    equal-contrib: true
    affiliation: 2
affiliations:
 - name: Computational Science Initiative, Brookhaven National Laboratory, Upton, New York 11973, USA
   index: 1
 - name: Center for Functional Nanomaterials, Brookhaven National Laboratory, Upton, New York 11973, USA
   index: 2
date: TODO
bibliography: paper.bib
# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Spectroscopy simulations are a critical tool for interpretation, the development of new theoretical understanding and fast screening of new molecules and materials. Systematically setting up input files for different simulation codes and multiple materials can be a time-consuming task with a relatively high barrier-to-entry, given the complexities and nuances of each individual simulation package. `Lightshow` solves this problem by providing a uniform abstraction for writing computational x-ray spectroscopy input files for multiple popular codes, including FEFF, VASP, OCEAN, EXCITING and Xspectra.

# Statement of need

<!-- `Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. -->

First-principles computer simulations of material and molecular properties are now a powerful tool at the forefront of bleeding-edge scientific research. On one hand, fundamental quantum mechanics equations can be solved numerically on high-performance computing systems, and the results interpreted at the microscopic level, something that is not possible otherwise. Thanks to their predictive nature, first-principles simulations provide fundamental understanding into the physical origins of various phenomena. On the other hand, they are critical to accelerating new materials design. Unlike expensive and time-consuming experiments, _in silico_ materials design frameworks can quickly screen the most promising candidates for target applications by running high-throughput calculations, allowing for a systematic down-sampling of the practically infinite chemical space. The emergence of high-performance computing hardware architecture combined with the development of efficient structure search algorithms continue to fuel the advance of first-principles simulations in domain science research.

Spectroscopy is an important experimental characterization technique that probes a sample based on the light-matter interaction. Different types of spectroscopy can be classified by the energy range they probe, such as X-ray, optical and infrared spectroscopy. In this work, we focus on one type: X-ray absorption spectroscopy (XAS), in which a deeply bound core level electron is excited to empty states, usually in the conduction bands. XAS is particularly useful because it is element-specific and very sensitive to local chemical environment of the absorbing sites, such as coordination number, charge state, and local symmetry. It has been widely used in condensed matter physics, geophysics, chemistry, materials science and biology for materials characterization. Recent instrument development at synchrotron light sources further improves the spatial, temporal and energy resolution of XAS, which opens new avenues in XAS research.

Despite the growing demand in first-principles XAS spectroscopy, carrying out practical calculations correctly is far from trivial, and requires a great deal of expertise in electronic structure theory, creating a formidable barrier for non-expert researchers. Most of the practical challenges boil down to the right choice of a set of input parameters, which depends on the level of theory, details of the implementation of the simulation software and the atomic structure of the system. A general purpose software package for generating XAS simulation input files for multiple codes does not exist. `Lightshow` is developed to fill this gap. It provides not only sets of default input parameters based on a careful multi-code XAS benchmark project (ref), but also exposes the entire suite of possible parameter choices for expert users to tune. Our goal is to provide an easy-to-use tool to the XAS community (for both newcomers and experts) for routine XAS simulation and analysis.

# Software description

![Caption for example figure.\label{fig:WorkflowDiagram}](figures/Lightshow_Workflow_Diagram.pdf)

We summarize the structure of `Lightshow` application programming interface (API) in Fig. \autoref{fig:WorkflowDiagram}. `Lightshow`'s core design philosophy is built around two principal objects, the `Database` class, and objects that inherit the `_BaseParameters` class. At a high-level, the `Database` class interfaces primarily with Pymatgen and the Materials Project `@Jain2013`, allowing the user to easily utilize the Pymatgen API and pull large numbers of materials structures quickly. Functionality is also available for instantiating a `Database` via loading e.g. `POSCAR`-style structure files from disk. Once a database has been created, code-specific simulation parameters inheriting the `_BaseParameters` base class interface with various methods in Pymatgen as well as in-house built software for systematically writing input files for multiple XAS simulation programs, including FEFF `@rehr2010parameter`, XSpectra `[@taillefumier2002x; @gougoussis2009first; bunuau2013projector]`, OCEAN `[@vinson2011bethe; @ocean-3]`, Exciting (ref) and VASP (ref). We highlight the `_tests` sub-directory, in which we maintain a suite of unit tests for individual functions and methods, as well as integration tests for the entire workflow.





<!-- 
# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:WorkflowDiagram}](figures/Lightshow_Workflow_Diagram.pdf)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }
 -->

# Acknowledgements

This material is based upon work supported by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, under Award Numbers FWP PS-030. The research used the theory and computation resources of the Center for Functional Nanomaterials, which is a U.S. DOE Office of Science Facility, and the Scientific Data and Computing Center, a component of the Computational Science Initiative, at Brookhaven National Laboratory under Contract No. DE-SC0012704.

# References

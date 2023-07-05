---
title: 'Lightshow: a Python package for generating computational x-ray absorption spectroscopy input files'
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
    orcid: 0000-0001-9869-9005
    equal-contrib: true
    corresponding: true
    affiliation: 2
  - name: Christian Vorwerk
    orcid: 0000-0002-2516-9553
    affiliation: 3
  - name: Benedikt Maurer
    orcid: 0000-0001-9152-7390
    affiliation: 4
  - name: Fabian Peschel
    orcid: 0000-0003-0619-6713
    affiliation: 4
  - name: Xiaohui Qu
    orcid: 0000-0001-5651-8405
    affiliation: 2
  - name: Eli Stavitski
    orcid: 0000-0002-3337-2930
    affiliation: 5
  - name: Claudia Draxl
    orcid: 0000-0003-3523-6657
    affiliation: 4
  - name: John Vinson
    orcid: 0000-0002-7619-7060
    affiliation: 6
  - name: Deyu Lu
    orcid: 0000-0003-4351-6085
    affiliation: 2
affiliations:
 - name: Computational Science Initiative, Brookhaven National Laboratory, Upton, New York 11973, United States
   index: 1
 - name: Center for Functional Nanomaterials, Brookhaven National Laboratory, Upton, New York 11973, United States
   index: 2
 - name: Pritzker School of Molecular Engineering, University of Chicago, Chicago, IL 60637, United States
   index: 3
 - name: Physics Department and IRIS Adlershof, Humboldt-Universität zu Berlin, D-12489 Berlin, Germany
   index: 4
 - name: National Synchrotron Light Source II, Brookhaven National Laboratory, Upton, New York 11973, United States
   index: 5
 - name: Material Measurement Laboratory, National Institute of Standards and Technology, Gaithersburg, Maryland 20899, United States
   index: 6
date: TODO
bibliography: paper.bib
# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

First-principles computational spectroscopy is a critical tool for interpreting experiment, 
performing structure refinement, and developing new physical understanding.
Systematically setting up input files for different simulation 
codes and a diverse class of materials is a challenging task with a very high 
barrier-to-entry, given the complexities and nuances of each individual simulation 
package. This task is non-trivial even for experts in the electronic structure field and nearly 
formidable for non-expert researchers.
`Lightshow` solves this problem by providing a uniform abstraction for 
writing computational x-ray spectroscopy input files for multiple popular codes, including 
FEFF, VASP, `OCEAN`, `exciting` and XSpectra. 
Its extendable framework will also allow the community to easily add new functions and 
to incorporate new simulation codes.

# Statement of need

First-principles simulations explore materials and molecular properties of a system by 
solving fundamental quantum mechanical equations numerically.
Thanks to their predictive nature, first-principles simulations provide fundamental understanding into the physical 
origins of various phenomena at the microscopic level, making them a powerful tool at the forefront
of a wide range of scientific research fields, including physics, chemistry, materials science and biology.
They are also critical to accelerating new materials 
design. In comparison to experiments, which can be expensive and time-consuming, _in silico_ materials 
design frameworks can quickly screen the most promising candidates for target 
applications by running high-throughput calculations, allowing for a systematic 
down-sampling of the intractably large chemical space. The emergence of
high-performance computing hardware architecture combined with the development 
of efficient structure search algorithms continues to fuel the advance of
first-principles simulations in materials design.

Spectroscopy is an important experimental characterization technique that
probes a sample based on the physics of light-matter interaction. Different types of
spectroscopy can be classified by the energy range they probe, such as X-ray,
ultraviolet–-visible, and infrared spectroscopy. `Lightshow` currently focuses on writing
the input files for one type: X-ray absorption spectroscopy (XAS), in which a
deeply bound core level electron is excited to empty states in the 
conduction bands. XAS is particularly useful because it is element-specific
and very sensitive to the local chemical environment of the absorbing sites,
such as coordination number, charge state, and local symmetry [@de2008core]. It has been
widely used in condensed matter physics, geophysics, chemistry, materials
science and biology for materials characterization. Recent instrument
development at synchrotron light sources further improves the spatial, temporal
and energy resolution of XAS, which opens new avenues in XAS research.

Despite the growing demands for first-principles XAS spectroscopy, carrying out
practical calculations correctly is far from trivial and requires a great deal
of expertise in electronic structure theory, creating a formidable barrier for
non-expert researchers. Most of the practical challenges boil down to the proper
choice of input parameters, which depends on the level of theory, details of
the implementation of the simulation software, and the atomic structure of the
system. A general purpose software package for generating XAS simulation input
files for multiple codes does not exist. `Lightshow` has been developed to
fill this gap. It provides not only sets of default input parameters based on
a careful multi-code XAS benchmark project [@xanesbench], but also exposes the
entire suite of possible parameter choices to expert users to tune. Our goal
is to provide an easy-to-use tool to the XAS community (for both newcomers
and experts) for XAS simulation and analysis. This tool will also help to improve
the data consistency and data reproducibility in the computational x-ray spectroscopy field, 
which are essential to data-driven applications.

# Brief software description

![Graphical representation of the organization of the `Lightshow` 
repository.\label{fig:WorkflowDiagram}](figures/Lightshow_Workflow_Diagram_2.pdf)

We summarize the structure of `Lightshow`'s application programming interface
(API) in \autoref{fig:WorkflowDiagram}. `Lightshow`'s core design
philosophy is built around two principal objects: the `Database` class and
the `_BaseParameters` class. At a high level, the `Database` class 
interfaces primarily with Pymatgen and the Materials Project [@jain2011high; @Jain2013], 
allowing the user to easily utilize Pymatgen and pull a large number of materials 
structures quickly. Functionality is also available for instantiating a 
`Database` via loading, e.g., `POSCAR`-style structure files from local 
storage. Once a database has been created, code-specific simulation parameters 
inheriting the `_BaseParameters` base class interface with various methods in 
Pymatgen as well as in-house built software for systematically writing input 
files for multiple XAS simulation programs, including FEFF [@rehr2010parameter; @kas2021advanced],
XSpectra [@taillefumier2002x; @gougoussis2009first; @bunuau2013projector],
`OCEAN` [@vinson2011bethe; @ocean-3], 
`exciting` [@exciting] and 
VASP [@vasp-xas].
We highlight the `_tests` directory, in which we maintain a suite of unit tests
for individual functions and methods, as well as integration tests for the
entire workflow. `Lightshow` is fully documented, and contains a simple example
notebook for users to get started. Finally, we also note that `Lightshow` is
designed to be part of larger workflows (including perhaps systematic comparison
to experiment) in which users wish to abstract away
the complicated and tedious task of generating input files. Future work on the
code will include designing modules for pre- and post-processing, allowing
for more seamless integration into said workflows.


# Acknowledgements

This work is partially supported by the U.S. Department of
Energy, Office of Science, Office of Basic Energy Sciences, under Award Numbers
FWP PS-030. The research used the theory and computation resources of the
Center for Functional Nanomaterials, which is a U.S. DOE Office of Science
Facility, and the Scientific Data and Computing Center, a component of the
Computational Science Initiative, at Brookhaven National Laboratory under
Contract No. DE-SC0012704. This work received partial funding by the German
Research Foundation (DFG) through the CRC 1404 (FONDA), Projektnummer
414984028, and the NFDI consortium FAIRmat – project 460197019.

C.V. acknowledges support by the Department of Energy, Basic Energy Sciences,
Materials Science and Engineering Division, through the Midwest Integrated
Center for Computational Materials (MICCoM).

Certain software is identified in this paper for clarity. Such identification
is not intended to imply recommendation or endorsement by NIST, nor is it
intended to imply that the materials or equipment identified are necessarily
the best available for the purpose.

# References

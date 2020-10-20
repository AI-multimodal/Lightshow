# J Vinson 2020
# C Vorwerk 2020
"""
Make the exciting input file and the photon files
"""

from ase.atoms import Atoms
from ase.io import write
import spglib
import pathlib
from os import environ as env
import sys

import pymatgen as pm
from pymatgen.io.exciting import ExcitingInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import numpy as np
#from photonSym import photonSymm
import json




def makeExciting( mpid, atoms: Atoms, params: dict ):
    
    # read default grst parameters from JSON file
    with open ("EXCITING/exciting.json", 'r') as fd:
        excitingJSON = json.load(fd)
    # read default xas parameters from JSON file
    with open ('EXCITING/exciting_xas.json', 'r') as fd:
        excitingxasJSON= json.load(fd)
    
    # generate pymatgen structure from ASE atoms
    struct = pm.io.ase.AseAtomsAdaptor.get_structure(atoms)
    
    # generate ExcitingInput object from pymatgen structure
    excitinginput = ExcitingInput(struct)
   
    # generate reduced and symmetrized structure
    finder = SpacegroupAnalyzer(struct)
    struct_red = finder.get_primitive_standard_structure(international_monoclinic = True)
    struct_symm = finder.get_symmetrized_structure()

    # find inequivalent sites
    ineq_atoms=[]
    for atoms in struct_symm.equivalent_indices:
        if struct_symm.sites[atoms[0]].species_string == params['species']:
            ineq_atoms.append(atoms[0])
   
    # set absorbing  edge
    excitingxasJSON['xs']['BSE']['xasedge'] = params['edge']
    
    # determine xasspecies parameter
    i=0
    for species in (sorted(struct_red.types_of_species, key=lambda el: el.X)):
        i=i+1
        if species.symbol == params['species']:
            excitingxasJSON['xs']['BSE']['xasspecies']=str(i)
    


    # generate filepath
    folder = pathlib.Path(env['PWD']) / mpid / "EXCITING"
    folder.mkdir(parents=True, exist_ok=True)
    filepath_grst=str(folder / 'input_grst.xml')
    
    # write input file for electronic-structure calculation
    try:
        excitinginput.write_file('primitive', filepath_grst, 
                                 bandstr=True, **excitingJSON)
    except:
        raise Exception("FAILED while trying to write input.xml")

    # write input file for XAS calculation
    for i in ineq_atoms:
    
        filepath_xas=str(folder / 'input_xas')+str(i+1)+'.xml'
        excitingxasJSON['xs']['BSE']['xasatom']=str(i+1)

        try:
            excitinginput.write_file('primitive', filepath_xas, 
                                     bandstr=False, **excitingxasJSON)
        except:
            raise Exception("FAILED while trying to write input.xml")

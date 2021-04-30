# J Vinson 2020
# C Vorwerk 2020
"""
Make the exciting input file and the photon files
"""

from ase.atoms import Atoms
from ase.io import write
import spglib
import pathlib
import os
from os import environ as env
import sys

import pymatgen as pm
from pymatgen.io.exciting import ExcitingInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import numpy as np
import json

import xanes_bench



def makeExcitingXAS( mpid, atoms: Atoms, params: dict ):
    
    # read default xas parameters from JSON file
    param_fn = os.path.join(os.path.dirname(xanes_bench.__file__), "EXCITING",
            "exciting_xas.json")
    with open (param_fn, 'r') as fd:
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
    folder = pathlib.Path(env['PWD']) / "data" / "mp_structures" / mpid / "EXCITING"
    folder.mkdir(parents=True, exist_ok=True)
    
    # write input file for XAS calculation
    for i in ineq_atoms:
    
        filepath_xas=str(folder / 'input_xas')+str(i+1)+'.xml'
        excitingxasJSON['xs']['BSE']['xasatom']=str(i+1)

        try:
            excitinginput.write_file('primitive', filepath_xas, 
                                     bandstr=False, **excitingxasJSON)
        except:
            raise Exception("FAILED while trying to write input.xml")

def makeExcitingGRST( mpid, atoms: Atoms, kpoints: list, nempty: int, filepath):
    
    # read default grst parameters from JSON file
    param_fn = os.path.join(os.path.dirname(xanes_bench.__file__), "EXCITING",
            "exciting.json")
    with open (param_fn, 'r') as fd:
        excitingJSON = json.load(fd)
    # set kpoints
    excitingJSON['groundstate']['ngridk']=" ".join([str(entry) for entry in kpoints])
    excitingJSON['groundstate']['nempty']=str(nempty)
    excitingJSON['structure']={}
    excitingJSON['structure']['autormt']="True"
    # generate pymatgen structure from ASE atoms
    struct = pm.io.ase.AseAtomsAdaptor.get_structure(atoms)
    
    # generate ExcitingInput object from pymatgen structure
    excitinginput = ExcitingInput(struct)
   
    # generate reduced and symmetrized structure
    finder = SpacegroupAnalyzer(struct)
    struct_red = finder.get_primitive_standard_structure(international_monoclinic = True)
    struct_symm = finder.get_symmetrized_structure()

    # generate filepath
    filepath_grst=str(filepath / 'input_grst.xml')
    
    # write input file for electronic-structure calculation
    try:
        excitinginput.write_file('primitive', filepath_grst, 
                                 bandstr=True, **excitingJSON)
    except:
        raise Exception("FAILED while trying to write input.xml")

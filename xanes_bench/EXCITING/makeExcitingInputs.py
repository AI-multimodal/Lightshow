# F Meng 2022
# J Vinson 2020
# C Vorwerk 2020
"""
Make the exciting input file and the photon files
"""

from pymatgen.core import Structure
import spglib
from pathlib import Path
import sys

import pymatgen as pm
from pymatgen.io.exciting import ExcitingInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

import numpy as np
import json

import xanes_bench

module_path = Path(xanes_bench.EXCITING.__path__[0])

def makeExcitingXAS( mpid, structure: Structure, params: dict ):
    ''' construct Exciting input files

        Parameters  
        ----------
        mpid : string, mandatory
            materials id from Materials Project
        structure : pymatgen.core.Structure, mandatory
            structure file
        params : dict, mandatory
            TODO

        Returns
        -------
        None
        * save input files to the corresponding directories
    '''
    # read default xas parameters from JSON file
    param_fn = module_path / "exciting_xas.json"
    with open (param_fn, 'r') as fd:
        excitingxasJSON= json.load(fd)
    
    # generate ExcitingInput object from pymatgen structure
    excitinginput = ExcitingInput(structure)
   
    # generate reduced and symmetrized structure
    finder = SpacegroupAnalyzer(structure)
    struct_red = finder.get_primitive_standard_structure(international_monoclinic = True)
    struct_symm = finder.get_symmetrized_structure()

    # find inequivalent sites
    ineq_atoms=[]
    for atoms in struct_symm.equivalent_indices:
        if struct_symm.sites[atoms[0]].species_string == params['species']:
            ineq_atoms.append(atoms[0])
   
    # set absorbing  edge
    excitingxasJSON['xs']['BSE']['xasedge'] = params['edge']

    # set BSE bands
    if params['conductionBands'] is not None:
        excitingxasJSON['xs']['BSE']['nstlxas'] = f"1 {params['conductionBands']}"

    # At the moment setting all k-point grids to be the same
    if params['scf.kpoints'] is not None:
        excitingxasJSON['groundstate']['ngridk'] = \
            f"{params['scf.kpoints'][0]} {params['scf.kpoints'][1]} {params['scf.kpoints'][2]}"
        excitingxasJSON['xs']['ngridk'] = \
            f"{params['scf.kpoints'][0]} {params['scf.kpoints'][1]} {params['scf.kpoints'][2]}"
        excitingxasJSON['xs']['ngridq'] = \
            f"{params['scf.kpoints'][0]} {params['scf.kpoints'][1]} {params['scf.kpoints'][2]}"
    
    # determine xasspecies parameter
    i=0
    for species in (sorted(struct_red.types_of_species, key=lambda el: el.X)):
        i=i+1
        if species.symbol == params['species']:
            excitingxasJSON['xs']['BSE']['xasspecies']=str(i)
    


    # generate filepath
    folder = Path.cwd() / "data" / "mp_structures" / mpid / "EXCITING" / \
             Path(f"Spectra-{params['scf.kpoints'][0]}-{params['scf.kpoints'][1]}-{params['scf.kpoints'][2]}")
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

def makeExcitingGRST( mpid, struct: Structure, kpoints: list, nempty: int, filepath):
    '''TODO
    '''
    # read default grst parameters from JSON file
    param_fn = module_path / "exciting.json"
    with open (param_fn, 'r') as fd:
        excitingJSON = json.load(fd)
    # set kpoints
    excitingJSON['groundstate']['ngridk']=" ".join([str(entry) for entry in kpoints])
    excitingJSON['groundstate']['nempty']=str(nempty)
    excitingJSON['structure']={}
    excitingJSON['structure']['autormt']="true"
    
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

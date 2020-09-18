# J Vinson 2020
"""
Make the ocean input file and the photon files
Later add in something to write out the pseudopotentials
"""

from ase.atoms import Atoms
from ase.io import write
import spglib
import pathlib
from os import environ as env
import sys
import numpy as np
#from photonSym import photonSymm
import json




def makeExciting( mpid, atoms: Atoms, params: dict ):


#    with open ("EXCITING/exciting.json", 'r') as fd:
#        excitingJSON = json.load(fd)

    
    us = {}
    ph = []
    photonSymm( atoms, us, ph, params['photonOrder'])


#    Right now we'll just calculate every atom for ocean instead of just the symmetry ones
#    oceanJSON['edges'] = ""
#    for site,  in us

#    symm= spglib.get_symmetry((atoms.get_cell(),
#                             atoms.get_scaled_positions(),
#                             atoms.get_atomic_numbers()),
#                             symprec=0.1, angle_tolerance=15)
#
#    equiv = symm['equivalent_atoms']

    symbols = atoms.get_chemical_symbols()

#    oceanJSON['toldfe'] = params['defaultConvPerAtom'] * len( symbols )

    
    folder = pathlib.Path(env['PWD']) / mpid / "EXCITING"
    folder.mkdir(parents=True, exist_ok=True)

    try:
        write(str(folder / "input.xml"), atoms, format='exciting' )
    except:
        raise Exception("FAILED while trying to write input.xml")


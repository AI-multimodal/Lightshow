# J Vinson 2020
"""
Make the ocean input file and the photon files
Later add in something to write out the pseudopotentials
"""

#from ase.atoms import Atoms
#from ase.io import write
from pymatgen.core import Structure
#from pymatgen.io.pwscf import PWInput

import spglib
import pathlib
from os import environ as env
import sys
import numpy as np
from xanes_bench.photonSym import photonSymm
import json
from xanes_bench.OCEAN.fakeASE import write_ocean_in
import xanes_bench.OCEAN
import os
from xanes_bench.General.kden import printKgrid

module_path = os.path.dirname(xanes_bench.OCEAN.__file__)

def makeOcean( mpid, structure: Structure, params: dict ):

    xs_fn = os.path.join(module_path, 'ocean.json')
    with open (xs_fn, 'r') as fd:
        oceanJSON = json.load(fd)
#    with open ("OCEAN/ocean.json", 'r') as fd:
#        oceanJSON = json.load(fd)

    if params['diemac'] is not None:
        oceanJSON['diemac'] = params['diemac']

    if params['conductionBands'] is not None:
        oceanJSON['nbands'] = -1 * params['conductionBands']

    if params['scf.kpoints'] is not None:
        oceanJSON['ngkpt'] = "{:d} {:d} {:d}".format( params['scf.kpoints'][0], params['scf.kpoints'][1], params['scf.kpoints'][2] )
    
    us = {}
    ph = []
    photonSymm( structure, us, ph, params['photonOrder'])


#    Right now we'll just calculate every atom for ocean instead of just the symmetry ones
#    oceanJSON['edges'] = ""
#    for site,  in us

#    symm= spglib.get_symmetry((atoms.get_cell(),
#                             atoms.get_scaled_positions(),
#                             atoms.get_atomic_numbers()),
#                             symprec=0.1, angle_tolerance=15)
#
#    equiv = symm['equivalent_atoms']

    #symbols = atoms.get_chemical_symbols()
    symbols = [str(i).split()[-1] for i in structure.species]
    oceanJSON['toldfe'] = params['defaultConvPerAtom'] * len( symbols )
    
    folder = pathlib.Path(env["PWD"]) / "data" / "mp_structures" / mpid / "OCEAN" / "Spectra"
    folder.mkdir(parents=True, exist_ok=True)
    printKgrid( structure, folder )
    print("length structure", len(structure))
    try:
        write_ocean_in(str(folder / "ocean.in"), structure, input_data=oceanJSON )
    except:
        raise Exception("FAILED while trying to write ocean.in")

    totalweight = 0
    for photon in ph:
        totalweight += photon["dipole"][3]

    photonCount = 0
    for photon in ph:
        photonCount += 1
        dir1 = photon["dipole"][0:3]
        dir2 = dir1
        weight = photon["dipole"][3] / totalweight
        mode = "dipole"

        with open( folder / ("photon%d" % (photonCount )), "w" ) as f:
            f.write( mode +"\n" )
            f.write( "cartesian %f %f %f \n" % (dir1[0], dir1[1], dir1[2] ) )
            f.write( "end\n" )
            f.write( "cartesian %f %f %f \n" % (dir2[0], dir2[1], dir2[2] ) )
            f.write( "end\n" )
            f.write( "4966\n" )  ### NEED TO FIX THIS (probably by moving it to a lookup table inside OCEAN)
            f.write( str(weight ) + "\n" )
            f.close

    # New total weight for the quadrupole terms
    totalweight = 0
    for photon in ph:
        if 'quad' in photon:
            for quad in photon["quad"]:
                totalweight += quad[3]

    for photon in ph:
        if 'quad' in photon:
            for quad in photon["quad"]:
                photonCount += 1 
                dir1 = photon["dipole"][0:3]
                dir2 = quad[0:3]
                weight = quad[3] / totalweight
                mode = "quad"

                with open( folder / ("photon%d" % (photonCount )), "w" ) as f:
                    f.write( mode +"\n" )
                    f.write( "cartesian %f %f %f \n" % (dir1[0], dir1[1], dir1[2] ) )
                    f.write( "end\n" )
                    f.write( "cartesian %f %f %f \n" % (dir2[0], dir2[1], dir2[2] ) )
                    f.write( "end\n" )
                    f.write( "4966\n" )  ### NEED TO FIX THIS (probably by moving it to a lookup table inside OCEAN)
                    f.write( str(weight ) + "\n" )
                    f.close

      

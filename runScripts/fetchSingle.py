# coding: latin-1
import numpy as np

from ase.atoms import Atoms
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor as ase

import sys
from math import exp

from Xspectra.makeXspectraInputs import makeXspectra
from OCEAN.makeOceanInputs import makeOcean
from EXCITING.makeExcitingInputs import makeExciting

# The script takes a single, positive integer to grab a system from materials project
if len(sys.argv) < 2 :
    print( "Requires MP number" )
    exit()
else :
    print( str(sys.argv) )
    mpid = 'mp-' + sys.argv[1]

# Your hashe materials project key needs to be in a file called mp.key
with open('mp.key', 'r' ) as f:
    mpkey = f.read()
    mpkey = mpkey.strip()

# MP api handler
#print( str(mpkey) )
mp = MPRester( str(mpkey) )


# These params should be off-loaded into a json file
params = dict(defaultConvPerAtom=1E-10, photonOrder=6)

# Build the unit cell within ASE by querying the material project api
unitC = ase.get_atoms(mp.get_structure_by_material_id(mpid, conventional_unit_cell=False))

#
# Grab some additional data from materials project (if it exists )
data = mp.query(criteria={"task_id": mpid}, properties=["diel","band_gap"])
#print( data[0] )

## Update OCEAN dielectric constant with calculated value or band_gap inverse-like
if  data[0]['diel'] is not None: 
    if data[0]['diel']['poly_electronic'] is not None:
        print( data[0]['diel']['poly_electronic'] )
        params['diemac'] = data[0]['diel']['poly_electronic']
        if data[0]['band_gap'] is not None:
            print( data[0]['band_gap'] )
            print( exp( 3.5/data[0]['band_gap'] ) )
elif data[0]['band_gap'] is not None:
    print( data[0]['band_gap'] )
    params['diemac'] = exp( 3.5/data[0]['band_gap'] )
    print(params['diemac'])




makeXspectra( mpid, unitC, params )

makeOcean( mpid, unitC, params )

makeExciting( mpid, unitC, params )


exit()

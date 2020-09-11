# coding: latin-1
import numpy as np
import pathlib
import spglib

from ase.atoms import Atoms
from ase.io import write
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor as ase

import sys
from OCEAN.fakeASE import write_ocean_in
import json

from pprint import pprint
from math import exp
import time
from pymatgen.ext.matproj import TaskType

from photonSym import photonSymm
from Xspectra.makeXspectraInputs import makeXspectra
from OCEAN.makeOceanInputs import makeOcean
from EXCITING.makeExcitingInputs import makeExciting


def main():

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

    st = mp.get_structure_by_material_id(mpid, conventional_unit_cell=False)
    st_dict = st.as_dict().copy()
    st_dict["download_at"] = time.ctime()
    with open(f"../data/mp_structures/{mpid}.json", 'w') as f:
        json.dump(st_dict, f, indent=4, sort_keys=True)
    unitC = ase.get_atoms(st)

    data = mp.query(criteria={"task_id": mpid}, properties=["diel","band_gap"])
    print( data[0] )

    params = dict(defaultConvPerAtom=1E-10, photonOrder=6)

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


if __name__ == '__main__':
    main()

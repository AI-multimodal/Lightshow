# coding: latin-1

import json
import sys
import time
from math import exp
import os

#from xanes_bench.EXCITING.makeExcitingInputs import makeExciting
# TODO: add implementation of photonSym to avoid import error
from xanes_bench.OCEAN.makeOceanInputs import makeOcean
from xanes_bench.Xspectra.makeXspectraInputs import makeXspectra
from xanes_bench.EXCITING.makeExcitingInputs import makeExcitingXAS

from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor as ase

import xanes_bench

# Return a guess at the number of conduction bands that a given unit-cell volume needs to 
# cover a given energy range (in Ryd)
#JTV
#TODO needs to be unified with groundState.py, both call the same function 
#     and the prefactor should always be the same
def getCondBands( volume, eRange):
    return round( 0.256 * volume * ( eRange**(3/2) ) )

def main():

    # The script takes a single, positive integer to grab a system from materials project
    if len(sys.argv) < 2 :
        print( "Requires MP number" )
        exit()
    else :
        print( str(sys.argv) )
        mpid = 'mp-' + sys.argv[1]

    # Your hashe materials project key needs to be in a file called mp.key
    mpkey_fn = os.path.join(os.path.dirname(xanes_bench.__file__), "mp.key")
    with open(mpkey_fn, 'r' ) as f:
        mpkey = f.read()
        mpkey = mpkey.strip()

    # MP api handler
    #print( str(mpkey) )
    mp = MPRester( str(mpkey) )

    st = mp.get_structure_by_material_id(mpid, conventional_unit_cell=False)
    st_dict = st.as_dict().copy()
    st_dict["download_at"] = time.ctime()
    st_dict["created_at"] = mp.get_doc(mpid)["created_at"]
    json_dir = "data"
    for spec_type in ["XS", "OCEAN"]:
        json_fn = f"{json_dir}/mp_structures/{mpid}/{spec_type}/Spectra/{mpid}.json"
        if not os.path.exists(os.path.dirname(json_fn)):
            os.makedirs(os.path.dirname(json_fn))
        with open(json_fn, 'w') as f:
            json.dump(st_dict, f, indent=4, sort_keys=True)
    unitC = ase.get_atoms(st)

    data = mp.query(criteria={"task_id": mpid}, properties=["diel","band_gap"])
    print( data[0] )

    cBands = getCondBands( unitC.get_volume(), 2.25 )
    params = dict(defaultConvPerAtom=1E-10, photonOrder=6, conductionBands=cBands)

    ## Update OCEAN dielectric constant with calculated value or band_gap inverse-like
    if  data[0]['diel'] is not None:
        if data[0]['diel']['poly_electronic'] is not None:
            print( data[0]['diel']['poly_electronic'] )
            params['diemac'] = data[0]['diel']['poly_electronic']
            if data[0]['band_gap'] is not None:
                print( data[0]['band_gap'] )
                if data[0]['band_gap'] > 0.000001:
                    print( exp( 3.5/data[0]['band_gap'] ) )
    elif data[0]['band_gap'] is not None:
        print( data[0]['band_gap'] )
        if data[0]['band_gap'] > 0.000001:
            params['diemac'] = exp( 3.5/data[0]['band_gap'] )
        else:
            params['diemac'] = 1000000
        print(params['diemac'])
    
    ## Add absorbing species and edge to parameters
    params['species']='Ti'
    params['edge']='K'



    makeXspectra( mpid, unitC, params )

    makeOcean( mpid, unitC, params )

    makeExcitingXAS( mpid, unitC, params )


if __name__ == '__main__':
    main()

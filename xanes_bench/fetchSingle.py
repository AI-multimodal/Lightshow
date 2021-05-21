# coding: latin-1

import json
import sys
import time
from math import exp, atan, pi
import os

#from xanes_bench.EXCITING.makeExcitingInputs import makeExciting
# TODO: add implementation of photonSym to avoid import error
from xanes_bench.OCEAN.makeOceanInputs import makeOcean
from xanes_bench.Xspectra.makeXspectraInputs import makeXspectra, makeXspectraConv
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

def dielectricGuess( gap ):
    if gap < 0.00001:
        return 1000000
    return 1.0/(2*atan(.164475*gap)/pi)

def main():

    # The script takes a single, positive integer to grab a system from materials project
    print(len(sys.argv))
    if len(sys.argv) < 2:
        print( "Requires MP number" )
        exit()
    else :
        print( str(sys.argv) )
        mpid = 'mp-' + sys.argv[1]
        if len(sys.argv) == 2:
            typecalc = "single"
            print("Type of Run: {}".format(typecalc))
        elif len(sys.argv) == 3: # Only use the second input paramater "single" or "convergence"
            if sys.argv[2] == "convergence" or sys.argv[2] == "single":
                typecalc = sys.argv[2]
                print("Type of Run: {}".format(typecalc))
            else:
                print("Input for run type not supported. \nSuppurted type of run: single or convergence")
                exit()

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
    for spec_type in ["XS", "OCEAN", "EXCITING"]:
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
    

    # Grab and parse k-point information
    # TODO error checking
    try:
        taskid = mp.query( criteria = {'task_id': mpid}, properties =
                ['blessed_tasks'])[0]['blessed_tasks']['GGA Static']
    except Exception as e:
        print(e)
        print( "Failed to get task id\nStopping\n")
        exit()
    try:
        data = mp.get_task_data( taskid, prop="kpoints" )
    except Exception as e:
        print(e)
        print( "Failed to get kpoints data\nStopping\n" )

    # If MP data set is incomplete, fail gracefully
    if 'kpoints' not in data[0]:
        print( "Failed to get kpoints from task id :", taskid )
        exit()

    # Returns a vasp kpoint object, so we need to convert to a dict
    kpointDict  = data[0]['kpoints'].as_dict()

    # We'll need to figure out the other types of grids, I thought I also saw Gamma
    #if kpointDict['generation_style'] != 'Monkhorst' :
    #    print( "Requires Monkhorst scheme" )
    #    exit()

    kpoints = kpointDict['kpoints'][0]
    koffset = kpointDict['usershift']


    # Make sure k-point grid is reasonable
    if kpoints[0]*kpoints[1]*kpoints[2] < 1 or kpoints[0]*kpoints[1]*kpoints[2] > 1000000 :
        print( "Bad k-point grid! ", kpoints )
        exit()

    params['scf.kpoints'] = kpoints


    ## Add absorbing species and edge to parameters
    params['species']='Ti'
    params['edge']='K'


    if typecalc == "single":

        makeXspectra( mpid, unitC, params )

        makeOcean( mpid, unitC, params )

        makeExcitingXAS( mpid, unitC, params )
    else:
        makeXspectraConv(mpid,unitC,params) # how to transfer kpoints?


if __name__ == '__main__':
    main()

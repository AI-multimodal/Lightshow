# coding: latin-1
  
import json
import sys
import time
import os, shutil

from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor as ase
import xanes_bench

from ase.atoms import Atoms
from ase.io import write

import pathlib
from os import environ as env

module_path = os.path.dirname(xanes_bench.__file__)

def main():

    # The script takes a single, positive integer to grab a system from materials project
    if len(sys.argv) < 2 :
        print( "Requires MP number" )
        exit()
    else :
        print( str(sys.argv) )
        mpid = 'mp-' + sys.argv[1]

    params = dict(defaultConvPerAtom=1E-10)

    # put in psp database handle here

    # Your hashed materials project key needs to be in a file called mp.key
    mpkey_fn = os.path.join(os.path.dirname(xanes_bench.__file__), "mp.key")
    with open(mpkey_fn, 'r' ) as f:
        mpkey = f.read()
        mpkey = mpkey.strip()

    # MP api handler
    mp = MPRester( str(mpkey) )

    st = mp.get_structure_by_material_id(mpid, conventional_unit_cell=False)
    st_dict = st.as_dict().copy()
    st_dict["download_at"] = time.ctime()
    st_dict["created_at"] = mp.get_doc(mpid)["created_at"]
    json_dir = f"../data"
    if not os.path.exists(json_dir):
        json_dir = "data"
    json_fn = f"{json_dir}/mp_structures/{mpid}.json"
    if not os.path.exists(os.path.dirname(json_fn)):
        os.makedirs(os.path.dirname(json_fn))
    with open(json_fn, 'w') as f:
        json.dump(st_dict, f, indent=4, sort_keys=True)
    unitC = ase.get_atoms(st)
    
    
    # Grab and parse k-point information
    # This will return a list of task IDs, but then we need to pick out the "correct" one
#    data = mp.get_data( mpid, prop="task_ids" )
#    print( data[0] )
    # We will want to change this up to search for a specific calculation type associated with mpid
    data = mp.get_task_data( mpid, prop="kpoints" )
    # Returns a vasp kpoint object, so we need to convert to a dict
    kpointDict  = data[0]['kpoints'].as_dict()
    # We'll need to figure out the other types of grids, I thought I also saw Gamma
    if kpointDict['generation_style'] != 'Monkhorst' :
        print( "Requires Monkhorst scheme" )
        exit()
    kpoints = kpointDict['kpoints'][0]
    koffset = kpointDict['usershift']
    print( koffset )
    # Might need to parse this, but it looks like most use Gamma-centered grids
    print( type( koffset ) )

    folder = pathlib.Path(env['PWD']) / mpid / "XS"
    folder.mkdir(parents=True, exist_ok=True)

    # defaults 
#    xs_fn = os.path.join(module_path, 'qe.json')
    qe_fn = 'qe.json'
    with open (qe_fn, 'r') as fd:
        qeJSON = json.load(fd)

    symbols = unitC.get_chemical_symbols()

    qeJSON['QE']['electrons']['conv_thr'] = params['defaultConvPerAtom'] * len( symbols )

#    sssp_fn = os.path.join(module_path, 'SSSP_precision.json')
    sssp_fn = 'SSSP_precision.json'
    psp = dict()
    with open (sssp_fn, 'r' ) as pspDatabaseFile:
        pspDatabase = json.load( pspDatabaseFile )
    minSymbols = set( symbols )
    for symbol in minSymbols:
        print( symbol )
        print( pspDatabase[ symbol ]['filename'] )
        psp[symbol] = pspDatabase[ symbol ]['filename']
        if qeJSON['QE']['system']['ecutwfc'] < pspDatabase[ symbol ]['cutoff']:
            qeJSON['QE']['system']['ecutwfc'] = pspDatabase[ symbol ]['cutoff']
        if 'rho_cutoff' in pspDatabase[ symbol ]:
            if 'ecutrho' not in qeJSON['QE']['system'] :
                 qeJSON['QE']['system']['ecutrho'] = pspDatabase[ symbol ]['rho_cutoff']
            if qeJSON['QE']['system']['ecutrho'] < pspDatabase[ symbol ]['rho_cutoff'] :
                qeJSON['QE']['system']['ecutrho'] = pspDatabase[ symbol ]['rho_cutoff']

        shutil.copy(
            os.path.join(module_path, "..", "data", "pseudopotential", "xspectral", "neutral",
                         pspDatabase[ symbol ]['filename']),
            str(folder / pspDatabase[symbol]['filename'])
        )

    shutil.copy(
        os.path.join(module_path, "..", "data", "pseudopotential", "xspectral", "orbital",
                     "Ti.wfc"),
        str(folder / "Ti.wfc")
    )

    try:
        write(str(folder / "qs.in"), unitC, format='espresso-in',
            input_data=qeJSON['QE'], pseudopotentials=psp, kpts=kpoints)
    except:
        print(qeJSON['QE'], unitC, psp)
        raise Exception("FAILED while trying to write qe.in")



if __name__ == '__main__':
    main()

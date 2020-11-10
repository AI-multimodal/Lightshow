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

import base64, bz2, hashlib

import re
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
import numpy as np

module_path = os.path.dirname(xanes_bench.__file__)

# Return a guess at the number of conduction bands that a given unit-cell volume needs to 
# cover a given energy range (in Ryd)
def getCondBands( volume, eRange):
    return round( 0.256 * volume * ( eRange**(3/2) ) )

#TODO add types for a bit of help
def writeQE( unitC, st, folder, qe_fn, pspName, params, conductionBands, kpoints ):

    with open (qe_fn, 'r') as fd:
        qeJSON = json.load(fd)

    symbols = unitC.get_chemical_symbols()

    qeJSON['QE']['electrons']['conv_thr'] = params['defaultConvPerAtom'] * len( symbols )
    qeJSON['QE']['control']['pseudo_dir'] = "../"

    psp_fn = os.path.join(module_path, "pseudos", "data", pspName + ".json" )
    with open (psp_fn, 'r' ) as pspDatabaseFile:
        pspDatabase = json.load( pspDatabaseFile )
    
    psp_fn = os.path.join(module_path, "pseudos", "data", pspName + "_pseudos.json")
    with open ( psp_fn, 'r' ) as pspDatabaseFile:
        pspFullData = json.load( pspDatabaseFile )


    psp = dict()
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

        if pspDatabase[ symbol ]['filename'] not in pspFullData:
            print( "Incomplete psp database" )
            exit()

        pspString = bz2.decompress(base64.b64decode( pspFullData[pspDatabase[ symbol ]['filename']] ))
        print( 'Expected hash:  ' + pspDatabase[symbol]['md5'] )
        print( 'Resultant hash: ' + hashlib.md5( pspString ).hexdigest() )


        fileName = os.path.join( folder, "..", pspDatabase[ symbol ]['filename'] )
        with open( fileName, 'w' ) as f:
            f.write( pspString.decode("utf-8") )


    nelectron = 0
    for symbol in symbols:
        nelectron += pspDatabase[ symbol ]["Z_val"]

    print( "N electron: ", nelectron)
#    qeJSON['QE']['system']['nbnd'] = round( nelectron/2 + conductionBands )

    # Write SCF input
    try:
        write(str(folder / "scf.in"), unitC, format='espresso-in',
            input_data=qeJSON['QE'], pseudopotentials=psp, kpts=kpoints)
    except:
        print(qeJSON['QE'], unitC, psp)
        raise Exception("FAILED while trying to write scf.in")



    # Write NSCF input for DOS (regular k-point mesh)
    qeJSON['QE']['system']['nbnd'] = round( nelectron/2 + conductionBands )
    qeJSON['QE']['control']['calculation'] = 'nscf'
    qeJSON['QE']['control']['tstress'] = False
    qeJSON['QE']['control']['tprnfor'] = False
    # This is new option as of 6.6, but should failsafe to default in earlier
    qeJSON['QE']['control']['disk_io'] = 'nowf'

    try:
        write(str(folder / "nscf.in"), unitC, format='espresso-in',
            input_data=qeJSON['QE'], pseudopotentials=psp, kpts=kpoints)
    except:
        print(qeJSON['QE'], unitC, psp)
        raise Exception("FAILED while trying to write nscf.in")


    # Now we want the band struture version
    #TODO revist number of bands?
    #TODO This is a hack to get around lack of k-path support
    ## We don't pass a kpoint spec which should give us "K_POINTS gamma"
    try:
        write(str(folder / "nscf_temp.in"), unitC, format='espresso-in',
            input_data=qeJSON['QE'], pseudopotentials=psp )
    except:
        print(qeJSON['QE'], unitC, psp)
        raise Exception("FAILED while trying to write nscf_temp.in")

    ## Now, slurp in entire file and get rid of K_POINT
    with open ( str(folder / "nscf_temp.in"), 'r' ) as f:
        NSCFtemp = re.sub('K_POINTS\s+gamma', '', f.read(), flags=re.IGNORECASE)

    ## Now do k-path
    qeJSON['QE']['control']['calculation'] = 'bands'
    ## Might have multiple k-point paths, best to break them into separate files

    # Might want to set the tolerances for this
    finder = SpacegroupAnalyzer(st)
    # This should be ok since MP is returning the primitive 
    new_struct = finder.get_primitive_standard_structure(international_monoclinic=False)

    kpath = HighSymmKpath(new_struct)
    # kpath has two parts
    # 'kpoints' gives the names of each high-symmetry point and its location
    # 'kpath' is a list of lists. There can be several paths that aren't connected

    #TODO This will need to be unified with Exciting
    ## Right now going for a target reciprocal space division, but I think 
    ## we might want to set an integer number to divide out instead
    targetKpointSpacing = 0.05

    bMatrix = new_struct.lattice.reciprocal_lattice.matrix

    # Loop over the separate paths
    for i in range(len(kpath.kpath['path'])):
        KString = "K_POINTS crystal_b\n%i\n" % len(kpath.kpath['path'][i])
        # Loop within a path
        symbol = kpath.kpath['path'][i][0]
        prevCoords = kpath.kpath['kpoints'][symbol]
        prevCart = np.dot( bMatrix, prevCoords )
        kpointCount = []
        totKpointCount = 0
        for j in range(1,len(kpath.kpath['path'][i])):
            symbol = kpath.kpath['path'][i][j]
            coords = kpath.kpath['kpoints'][symbol]
            cart = np.dot( bMatrix, coords )
            dist = np.linalg.norm(cart-prevCart)
            kpointCount.append( int( dist/targetKpointSpacing ) )
    #        prevCoords = coords
            prevCart = cart
            totKpointCount += int( dist/targetKpointSpacing )

        kpointCount.append( int(1) )
#        print( totKpointCount )
#        print( len(kpath.kpath['path'][i]) )

        for j in range(len(kpath.kpath['path'][i])):
            symbol = kpath.kpath['path'][i][j]
            coords = kpath.kpath['kpoints'][symbol]
#            print( "%16.12f %16.12f %16.12f %i" % (coords[0],coords[1],coords[2],kpointCount[j]) )
            KString += "%16.12f %16.12f %16.12f %i\n" % (coords[0],coords[1],coords[2],kpointCount[j])

#        print( "  " )

        with open ( str(folder / "nscf_band" ) + ".%i.in" % (i+1), 'w' ) as f:
            f.write( NSCFtemp )
            f.write( KString )




    
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
    #TODO gracefully report errors with connection, fetching structure

    st_dict = st.as_dict().copy()
    st_dict["download_at"] = time.ctime()
    st_dict["created_at"] = mp.get_doc(mpid)["created_at"]
    json_dir = "data"
    for spec_type in ["XS", "OCEAN"]:
        json_fn = f"{json_dir}/mp_structures/{mpid}/{spec_type}/groundState/{mpid}.json"
        if not os.path.exists(os.path.dirname(json_fn)):
            os.makedirs(os.path.dirname(json_fn))
        with open(json_fn, 'w') as f:
            json.dump(st_dict, f, indent=4, sort_keys=True)
    unitC = ase.get_atoms(st)
    
    conductionBands = getCondBands( unitC.get_volume(), 1.5 )
    print( "Conduction bands: ", conductionBands )
    
    # Grab and parse k-point information
    # TODO error checking
    taskid = mp.query( criteria = {'task_id': mpid}, properties =
            ['blessed_tasks'])[0]['blessed_tasks']['GGA Static']
    data = mp.get_task_data( taskid, prop="kpoints" )

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


    # Right now we are not checking for grid shifts. QE will just use a Gamma-centered grid


    
    # defaults, will be common for both "ocean" and "XS" as they are both (for now) using QE
    qe_fn = os.path.join(module_path, 'QE', 'qe.json')

    # subdir says where to put the input and psps 
    subdir = pathlib.Path(env['PWD'], "data", "mp_structures", mpid, "XS", "groundState")
    subdir.mkdir(parents=True, exist_ok=True)
    writeQE( unitC, st, subdir , qe_fn, 'SSSP_precision', params, conductionBands, kpoints )


    subdir = pathlib.Path(env['PWD'], "data", "mp_structures",mpid, "OCEAN", "groundState" )
    subdir.mkdir(parents=True, exist_ok=True)
    writeQE( unitC, st, subdir , qe_fn, 'PD_stringent', params, conductionBands, kpoints )

    #TODO should be able to add in calls to exciting io here


if __name__ == '__main__':
    main()

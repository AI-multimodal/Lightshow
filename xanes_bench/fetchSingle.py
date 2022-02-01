# coding: latin-1
# Fanchen Meng, 2022
# based on John Vinson's version
import json
import sys
import numpy as np
from math import exp
from pathlib import Path

from xanes_bench.OCEAN.makeOceanInputs import makeOcean
from xanes_bench.Xspectra.makeXspectraInputs import makeXspectra
from xanes_bench.EXCITING.makeExcitingInputs import makeExcitingXAS
from xanes_bench.General.kden import returnKpoint
from xanes_bench.utils import * # TODO
from xanes_bench.Xspectra.preprocessing import build_supercell
import xanes_bench

from pymatgen.io.ase import AseAtomsAdaptor as ase

def main():
    # The script takes a single, positive integer to grab a system from materials project
    # determine the type of calculation
    if len(sys.argv) < 2:
        print( "Requires MP number" )
        exit()
    else:
        mpid = sys.argv[1]
        # determine type of calculations
        # default : single, can also be set explicitly 
        # TODO : others? how to do convergence?
        if len(sys.argv) == 2:
            typecalc = "single"
            print(f"Type of Run: {typecalc}")
        else:
            if sys.argv[2] == "single":
                typecalc = sys.argv[2]
                print(f"Type of Run: {typecalc}")
            else:
                print("Input for run type not supported. \nSuppurted type of run: single") 
                exit()

    # get structure from Materials Project 
    mpr = setMPR()
    st, st_dict = get_structure(mpid)
    if typecalc == "single":
        json_dir = "data"
        json_fn = Path(f"{json_dir}/mp_structures/{mpid}/{mpid}.json") 
        if not Path.exists(json_fn.parent):
            Path.mkdir(json_fn.parent, parents=True, exist_ok=True)
        with open(json_fn, 'w') as f:
            json.dump(st_dict, f, indent=4, sort_keys=True)

    data = mpr.query(criteria={"task_id": mpid}, properties=["diel","band_gap"])
    cBands = getCondBands( st.lattice.volume, 2.25 ) 
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

    # get k-points using 45 Bohr threshold
    # valid for unitcell calculations, e.g. OCEAN, EXCITING
    kpoints = returnKpoint(st, 23.813) # 23.813 A = 45 Bohr
    koffset = [0.0, 0.0, 0.0]
    params['scf.kpoints'] = kpoints

    ## Add absorbing species and edge to parameters
    params['species']='Ti'
    params['edge']='K'


    if typecalc == "single":
        makeOcean( mpid, st, params )
        makeExcitingXAS( mpid, st, params )
        # build supercell for xspectra here
        st = build_supercell(ase.get_atoms(st))
        kpoints = returnKpoint(st, 23.813)
        params['scf.kpoints'] = kpoints
        makeXspectra( mpid, st, params )

if __name__ == '__main__':
    main()

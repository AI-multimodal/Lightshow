# F Meng 2022
# J Vinson 2020
"""
Make the ocean input file and the photon files
Later add in something to write out the pseudopotentials
"""

from pymatgen.core import Structure

import spglib
from pathlib import Path
import numpy as np
import json
from xanes_bench.photonSym import photonSymm
from xanes_bench.OCEAN.fakeASE import write_ocean_in
import xanes_bench.OCEAN
from xanes_bench.General.kden import printKgrid

module_path = Path(xanes_bench.OCEAN.__path__[0])

def makeOcean( mpid, structure: Structure, params: dict ):
    ''' construct ocean input files

        Parameters
        ----------
        mpid : str, mandatory
            materials id in Materials Project
        structure : pymatgen.core.Structure, mandatory
            structure data
        params : dict, mandatory
            TODO

        Returns
        -------
        None
        * save input files to corresponding directories
    '''
    xs_fn = module_path / 'ocean.json'
    with open (xs_fn, 'r') as fd:
        oceanJSON = json.load(fd)

    if params['diemac'] is not None:
        oceanJSON['diemac'] = params['diemac']

    if params['conductionBands'] is not None:
        oceanJSON['nbands'] = -1 * params['conductionBands']

    if params['scf.kpoints'] is not None:
        oceanJSON['ngkpt'] = \
                f"{params['scf.kpoints'][0]} {params['scf.kpoints'][1]} {params['scf.kpoints'][2]}"
    
    us = {}
    ph = []
    photonSymm( structure, us, ph, params['photonOrder'])


#    Right now we'll just calculate every atom for ocean instead of just the symmetry ones
    symbols = [str(i).split()[-1] for i in structure.species]
    oceanJSON['toldfe'] = params['defaultConvPerAtom'] * len( symbols )
    json_dir = params['json_dir']
    folder = Path.cwd() / json_dir / "mp_structures" / mpid / "OCEAN" / \
            Path(f"Spectra-{params['scf.kpoints'][0]}-{params['scf.kpoints'][1]}-{params['scf.kpoints'][2]}")
    folder.mkdir(parents=True, exist_ok=True)
    printKgrid( structure, folder )
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

      

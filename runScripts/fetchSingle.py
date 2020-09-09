# coding: latin-1
import numpy as np
import pathlib
import spglib

from ase.atoms import Atoms
from ase.io import write
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor as ase
from os import environ as env
import sys
from OCEAN.fakeASE import write_ocean_in
import json

from pprint import pprint
from math import exp

from photonSym import photonSymm
from Xspectra.makeXspectraInputs import makeXspectra
from OCEAN.makeOceanInputs import makeOcean

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

psp = dict(Ti1='Ti.fch.upf')


###### ocean defaultes
ocean_input = { 'dft': 'qe', 'opf.program': 'hamann', 'para_prefix': 'mpirun -n 8',
                'pp_database': 'ONCVPSP-PBE-PDv0.4-stringent', 'edges': '-22 1 0',
                'ecut': '-1', 'diemac': '5', 'dft_energy_range': 50, 'screen_energy_range': 150, 
                'ecut.quality': 'high' }



unitC = ase.get_atoms(mp.get_structure_by_material_id(mpid, conventional_unit_cell=False))

data = mp.query(criteria={"task_id": mpid}, properties=["diel","band_gap"])
print( data[0] )

params = dict(defaultConvPerAtom=1E-10, photonOrder=6)

## Update OCEAN dielectric constant with calculated value or band_gap inverse-like
if  data[0]['diel'] is not None: 
    if data[0]['diel']['poly_electronic'] is not None:
        print( data[0]['diel']['poly_electronic'] )
#        ocean_input['diemac'] = data[0]['diel']['poly_electronic']
        params['diemac'] = data[0]['diel']['poly_electronic']
        if data[0]['band_gap'] is not None:
            print( data[0]['band_gap'] )
            print( exp( 3.5/data[0]['band_gap'] ) )
elif data[0]['band_gap'] is not None:
    print( data[0]['band_gap'] )
#    ocean_input['diemac'] = exp( 3.5/data[0]['band_gap'] )
#    print(ocean_input['diemac'])
    params['diemac'] = exp( 3.5/data[0]['band_gap'] )
    print(params['diemac'])




makeXspectra( mpid, unitC, params )

makeOcean( mpid, unitC, params )

###############
############## JTV !!!!
###############
exit()
###############
############## JTV !!!!
###############


## Set OCEAN conv threshold
conv_thr = defaultConvPerAtom * len( unitC.numbers )
ocean_input['toldfe'] = conv_thr

## Make supercell for QE/xspectra
atoms = smaller( unitC )

photonSymm( atoms )

## Get symm for supercell
symm= spglib.get_symmetry((atoms.get_cell(),
                           atoms.get_scaled_positions(),
                           atoms.get_atomic_numbers()),
                           symprec=0.1, angle_tolerance=15)

equiv = symm['equivalent_atoms']

symbols = atoms.get_chemical_symbols()

## Set QE conv thr per atom
conv_thr = defaultConvPerAtom * len( symbols )
electrons['conv_thr'] = conv_thr

## Look up psps from database, save names and update cut-off energies
pspDatabaseFile = open( 'SSSP_precision.json' )
pspDatabase = json.load( pspDatabaseFile )
minSymbols = set( symbols )
for symbol in minSymbols:
    print( symbol )
    print( pspDatabase[ symbol ]['filename'] )
    psp[symbol] = pspDatabase[ symbol ]['filename']
    if system['ecutwfc'] < pspDatabase[ symbol ]['cutoff']:
        system['ecutwfc'] = pspDatabase[ symbol ]['cutoff']
    if system['ecutrho'] < pspDatabase[ symbol ]['rho_cutoff']:
        system['ecutrho'] = pspDatabase[ symbol ]['rho_cutoff']
    

iabs_ = 0
iabs = []
found = set()

for symbol in symbols:
    if symbol == 'Ti':
        iabs.append(len(found) + 1)
    else:
        iabs.append(None)
    found.add(symbol)

prev = None

input = dict(
    control=control, system=system, electrons=electrons, ions=ions, cell=cell)

folder = pathlib.Path(env['PWD']) / mpid
folder.mkdir(parents=True, exist_ok=True)
try:
    write(str(folder / "qe.in"), atoms, format='espresso-in',
        input_data=input, pseudopotentials=psp, kpts=[1, 1, 1])
except:
    print(input, atoms, psp)
    raise Exception("FAILED while trying to write qe.in")

try:
    write_ocean_in(str(folder / "ocean.in"), unitC, input_data=ocean_input )
except:
    raise Exception("FAILED while trying to write ocean.in")

pprint( equiv )

for i, sym in enumerate(symbols):

    if i == equiv[i] and sym == 'Ti':

        if prev is not None:
            atoms[prev].tag = 0

        atoms[i].tag = 1
        prev = i

        subfolder = folder / str(i)
        subfolder.mkdir(parents=True, exist_ok=True)

        write(str(subfolder / "input.txt"), atoms, format='espresso-in',
            input_data=input, pseudopotentials=psp, kpts=[1, 1, 1])

        for dir in (1, 2, 3):
            dir1, dir2 = directions - set((dir, ))
            dirs = (dir, dir1, dir2)

            for mode in ("dipole", "quadrupole"):

                xanesfolder = subfolder / ("%s%d" % (mode, dir))

                if any([np.all(XY == rot) for rot in symm['rotations']]) and dir == 2:
                    xanesfolder.symlink_to(("%s1"%mode))
                    continue

                if any([np.all(XZ == rot) for rot in symm['rotations']]) and dir == 3:
                    xanesfolder.symlink_to(("%s1"%mode))
                    continue

                if any([np.all(YZ == rot) for rot in symm['rotations']]) and dir == 3:
                    xanesfolder.symlink_to(("%s2"%mode))
                    continue

                xanesfolder.mkdir(parents=True, exist_ok=True)

                with open(xanesfolder / "xanes.in", "w") as f:
                    f.write(xinput(mode, iabs[i], dirs,
                                   ortho(atoms.cell, dir - 1 ,
                                         (dir2 if dir==2 else dir1)-1)))

                with open(xanesfolder / "xanes_.in", "w") as f:
                    f.write(xinput(mode, iabs[i], dirs,
                                   ortho(atoms.cell, dir - 1 ,
                                         (dir2 if dir==2 else dir1)-1),
                                   plot=True))

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

### Sym ops for xspectra
XY = np.array(([0, 1, 0], [1, 0, 0], [0, 0, 1]))
XZ = np.array(([0, 0, 1], [0, 1, 0], [1, 0, 0]))
YZ = np.array(([1, 0, 0], [0, 0, 1], [0, 1, 0]))

directions = {1, 2, 3}

def smaller(atoms):
#    # starting from the primitive cell, give a supercell
#    # that has at least 9 Å in each direction; either using
#    # primitive or conventinal cell as building blocks;
#    # returns which ever is smaller
    prim =  atoms * ((9 / np.linalg.norm(atoms.cell, axis=1)).astype(int) + 1)
    lat, pos, Z = spglib.standardize_cell((atoms.get_cell(),
                                           atoms.get_scaled_positions(),
                                           atoms.get_atomic_numbers()))
    conv = Atoms(Z, cell=lat, positions=pos@lat, pbc=True)
    conv = conv * ((9 / np.linalg.norm(conv.cell, axis=1)).astype(int) + 1)
    return conv if len(conv) <= len(prim) else prim

# For both OCEAN and QE/xspectra we want to have a uniform level of convergence (per atom?)
defaultConvPerAtom = 1E-10

##### These dictionaries are for building QE/xspectra inputs
# input file description
#psp = dict(Ti='Ti.upf',
#           Ti1='Ti.fch.upf',
#           O='O.upf')
psp = dict(Ti1='Ti.fch.upf')

control = dict(restart_mode='from_scratch',
               wf_collect=True)

#set ecutwfc and ecutrho to match the fch psp
system = dict(
    ecutwfc=40,
    ecutrho=320,
    occupations='smearing',
    smearing='mv',
    degauss=0.002,
    nspin=1)

electrons = dict(mixing_beta=0.4, conv_thr=1E-8)

ions = dict()
cell = dict()
###### end QE/xspectra defaults

# Move lower to update ecut
#input = dict(
#    control=control, system=system, electrons=electrons, ions=ions, cell=cell)

###### ocean defaultes
ocean_input = { 'dft': 'qe', 'opf.program': 'hamann', 'para_prefix': 'mpirun -n 8',
                'pp_database': 'ONCVPSP-PBE-PDv0.4-stringent', 'edges': '-22 1 0',
                'ecut': '-1', 'diemac': '5', 'dft_energy_range': 50, 'screen_energy_range': 150, 
                'ecut.quality': 'high' }


def ortho(lat, i, j):
    o = np.linalg.solve(lat.T, np.cross(lat[i], lat[j]))
    return tuple(o / np.linalg.norm(o))

def xinput(mode, iabs, dirs, xkvec, plot=False):

    inp =  ["&input_xspectra",
            "    calculation = 'xanes_%s'" % mode,
            "    edge = 'K'",
            "    prefix = 'pwscf'",
            "    outdir = '../'",
            "    xniter = 5000",
            "    xiabs = %d" % iabs,
            "    xerror = 0.01",
            "    wf_collect = .true.",
            "    xcheck_conv = 200",
            "    xepsilon(%d) = 1.0" % dirs[0],
            "    xepsilon(%d) = 0.0" % dirs[1],
            "    xepsilon(%d) = 0.0" % dirs[2]]

    if mode == "quadrupole":

        inp += ["    xkvec(1) = %.10f" % xkvec[0],
                "    xkvec(2) = %.10f" % xkvec[1],
                "    xkvec(3) = %.10f" % xkvec[2] ]

    if plot:
        inp += ["    xonly_plot = .true."]

    inp += ["/",
            "&plot",
            "    xnepoint = 400",
            "    xemin = -5.0",
            "    xemax = 60",
            "    terminator = .true.",
            "    cut_occ_states = .true."]
    if plot:
        inp += [
            "    xgamma = 0.05",
            "    gamma_mode = 'constant'",
            "/"]
    else:
        inp += [
            "    gamma_mode = 'variable'",
            "    gamma_energy(1) = 7",
            "    gamma_energy(2) = 23",
            "    gamma_value(1)  = 0.89",
            "    gamma_value(2)  = 2.1",
            "/"]

    inp += ["&pseudos",
            "    filecore = '../../../Ti.wfc'",
            "/",
            "&cut_occ",
            "    cut_desmooth = 0.3",
            "/",
            "2 2 2 0 0 0"]
    return '\n'.join(inp) + '\n'



unitC = ase.get_atoms(mp.get_structure_by_material_id(mpid, conventional_unit_cell=False))

data = mp.query(criteria={"task_id": mpid}, properties=["diel","band_gap"])
print( data[0] )


## Update OCEAN dielectric constant with calculated value or band_gap inverse-like
if  data[0]['diel'] is not None: 
    if data[0]['diel']['poly_electronic'] is not None:
        print( data[0]['diel']['poly_electronic'] )
        ocean_input['diemac'] = data[0]['diel']['poly_electronic']
        if data[0]['band_gap'] is not None:
            print( data[0]['band_gap'] )
            print( exp( 3.5/data[0]['band_gap'] ) )
elif data[0]['band_gap'] is not None:
    print( data[0]['band_gap'] )
    ocean_input['diemac'] = exp( 3.5/data[0]['band_gap'] )
    print(ocean_input['diemac'])

## Set OCEAN conv threshold
conv_thr = defaultConvPerAtom * len( unitC.numbers )
ocean_input['toldfe'] = conv_thr


us = {}
ph = ()
#photonSymm( unitC, us, ph)

params = dict(defaultConvPerAtom=1E-10, photonOrder=6)

makeXspectra( mpid, unitC, params )

exit()
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

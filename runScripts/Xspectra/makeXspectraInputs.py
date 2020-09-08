# John Vinson, 2020
# based on earlier work by Sencer
"""
Make all the xspectra/qe inputs
"""
from ase.atoms import Atoms
import ase.io.espresso
from ase.io import write
import spglib
import pathlib
from os import environ as env
import sys
import numpy as np
from photonSym import photonSymm
import json


def smaller(atoms: Atoms, Rmin=9.0):
#    # starting from the primitive cell, give a supercell
#    # that has at least 9 Å in each direction; either using
#    # primitive or conventinal cell as building blocks;
#    # returns which ever is smaller
    prim =  atoms * ((Rmin / np.linalg.norm(atoms.cell, axis=1)).astype(int) + 1)
    lat, pos, Z = spglib.standardize_cell((atoms.get_cell(),
                                           atoms.get_scaled_positions(),
                                           atoms.get_atomic_numbers()))
    conv = Atoms(Z, cell=lat, positions=pos@lat, pbc=True)
    conv = conv * ((Rmin / np.linalg.norm(conv.cell, axis=1)).astype(int) + 1)
    return conv if len(conv) <= len(prim) else prim


def ortho(lat, i, j):
    o = np.linalg.solve(lat.T, np.cross(lat[i], lat[j]))
    return tuple(o / np.linalg.norm(o))

def xinput(mode, iabs, dirs, xkvec, XSparams: dict, plot=False):

    inp =  ["&input_xspectra",
            "    calculation = 'xanes_%s'" % mode,
            "    edge = '" + XSparams['input_xspectra']['edge'] + "'",
            "    prefix = 'pwscf'",
            "    outdir = '../'",
            "    xniter = " + str(XSparams['input_xspectra']['xniter']),
            "    xiabs = %d" % iabs,
            "    xerror = " + str(XSparams['input_xspectra']['xerror']),
            "    wf_collect = .true.",
            "    xcheck_conv = " + str(XSparams['input_xspectra']['xcheck_conv']),
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
            "    xnepoint = " + str(XSparams['plot']['xnepoint']),
            "    xemin = " + str(XSparams['plot']['xemin']),
            "    xemax = " + str(XSparams['plot']['xemax']),
            "    terminator = " + XSparams['plot']['terminator'],
            "    cut_occ_states = " + XSparams['plot']['cut_occ_states'] ]
    if plot:
        inp += [
            "    xgamma = " + str(XSparams['plotTrue']['xgamma']),
            "    gamma_mode = '" + XSparams['plotTrue']['gamma_mode'] + "'",
            "/"]
    else:
        inp += [
            "    gamma_mode = '" + XSparams['plotFalse']['gamma_mode'] + "'",
            "    gamma_energy(1) = " + str( XSparams['plotFalse']['gamma_energy(1)']),
            "    gamma_energy(2) = " + str( XSparams['plotFalse']['gamma_energy(2)']),
            "    gamma_value(1)  = " + str( XSparams['plotFalse']['gamma_value(1)']),
            "    gamma_value(2)  = " + str( XSparams['plotFalse']['gamma_value(2)']),
            "/"]

    inp += ["&pseudos",
            "    filecore = '../../../Ti.wfc'",
            "/",
            "&cut_occ",
            "    cut_desmooth = " + str( XSparams['cut_occ']['cut_desmooth']),
            "/",
            XSparams['kpts']['kpts'] + " " + XSparams['kpts']['shift'] ]
    return '\n'.join(inp) + '\n'



def makeXspectra( mpid, unitCell: Atoms, params: dict ):
    #######
    psp = dict(Ti1='Ti.fch.upf')
    symTarg = 'Ti'
    ####

    with open ("Xspectra/xspectra.json", 'r') as fd:
        xsJSON = json.load(fd)

    atoms = smaller( unitCell )

    us = {}
    ph = []
    photonSymm( atoms, us, ph, params['photonOrder'])

    print( ph )

    symm= spglib.get_symmetry((atoms.get_cell(),
                             atoms.get_scaled_positions(),
                             atoms.get_atomic_numbers()),
                             symprec=0.1, angle_tolerance=15)

    equiv = symm['equivalent_atoms']

    symbols = atoms.get_chemical_symbols()


    xsJSON['QE']['electrons']['conv_thr'] = params['defaultConvPerAtom'] * len( symbols )

    with open ('Xspectra/SSSP_precision.json', 'r' ) as pspDatabaseFile:
        pspDatabase = json.load( pspDatabaseFile )
    minSymbols = set( symbols )
    for symbol in minSymbols:
        print( symbol )
        print( pspDatabase[ symbol ]['filename'] )
        psp[symbol] = pspDatabase[ symbol ]['filename']
        if xsJSON['QE']['system']['ecutwfc'] < pspDatabase[ symbol ]['cutoff']:
            xsJSON['QE']['system']['ecutwfc'] = pspDatabase[ symbol ]['cutoff']
        if xsJSON['QE']['system']['ecutrho'] < pspDatabase[ symbol ]['rho_cutoff']:
            xsJSON['QE']['system']['ecutrho'] = pspDatabase[ symbol ]['rho_cutoff']



#    input = xsJSON['QE']
    folder = pathlib.Path(env['PWD']) / mpid / "XS"
    folder.mkdir(parents=True, exist_ok=True)
    try:
        write(str(folder / "qe.in"), atoms, format='espresso-in',
            input_data=xsJSON['QE'], pseudopotentials=psp, kpts=[1, 1, 1])
    except:
        print(xsJSON['QE'], atoms, psp)
        raise Exception("FAILED while trying to write qe.in")

    iabs_ = 0
    iabs = []
    found = set()

    for symbol in symbols:
        if symbol == symTarg:
            iabs.append(len(found) + 1)
        else:
            iabs.append(None)
        found.add(symbol)

    prev = None
    XY = np.array(([0, 1, 0], [1, 0, 0], [0, 0, 1]))
    XZ = np.array(([0, 0, 1], [0, 1, 0], [1, 0, 0]))
    YZ = np.array(([1, 0, 0], [0, 0, 1], [0, 1, 0]))

    directions = {1, 2, 3}


    for i, sym in enumerate(symbols):

      if i == equiv[i] and sym == symTarg:

          if prev is not None:
              atoms[prev].tag = 0

          atoms[i].tag = 1
          prev = i

          subfolder = folder / str(i)
          subfolder.mkdir(parents=True, exist_ok=True)

          write(str(subfolder / "input.txt"), atoms, format='espresso-in',
              input_data=xsJSON['QE'], pseudopotentials=psp, kpts=[1, 1, 1])

          # OCEAN photon labeling is continuous, so we will do that here too
          #  not sure that we will actually want dipole-only spectra(?)
          totalweight = 0
          for photon in ph:
              totalweight += photon["dipole"][3]
#          totalweight *= Somethingforthesymmetry
          
          photonCount = 0
          for photon in ph:
#              print( dipole )
              photonCount += 1
              dir1 = photon["dipole"][0:3]
              dir2 = dir1
              weight = photon["dipole"][3] * us[i] / totalweight
              mode = "dipole"
              xanesfolder = subfolder / ("%s%d" % (mode, photonCount))
              xanesfolder.mkdir(parents=True, exist_ok=True)
              # fo
              with open(xanesfolder / "xanes.in", "w") as f:
                      f.write(xinput(mode, iabs[i], dir1,
                                           dir2, xsJSON['XS']))

              with open(xanesfolder / "xanes_.in", "w") as f:
                  f.write(xinput(mode, iabs[i], dir1,
                                 dir2,  xsJSON['XS'], plot=True)) 

              with open(xanesfolder / "weight.txt", "w") as f:
                  f.write( str(weight) + "\n" )


          # New total weight for the quadrupole terms
          totalweight = 0
          for photon in ph:
              for quad in photon["quad"]:
                  totalweight += quad[3]

          for photon in ph:
              for quad in photon["quad"]:
                  photonCount += 1
                  dir1 = photon["dipole"][0:3]
                  dir2 = quad[0:3]
                  weight = quad[3] * us[i] / totalweight
                  mode = "quadrupole"
                  xanesfolder = subfolder / ("%s%d" % (mode, photonCount))
                  xanesfolder.mkdir(parents=True, exist_ok=True)
                  # fo
                  with open(xanesfolder / "xanes.in", "w") as f:
                          f.write(xinput(mode, iabs[i], dir1,
                                               dir2, xsJSON['XS']))

                  with open(xanesfolder / "xanes_.in", "w") as f:
                      f.write(xinput(mode, iabs[i], dir1,
                                     dir2,  xsJSON['XS'], plot=True))

                  with open(xanesfolder / "weight.txt", "w") as f:
                      f.write( str(weight) + "\n" )

"""              
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
                                     ortho(atoms.cell, dir - 1, 
                                           (dir2 if dir==2 else dir1)-1), xsJSON['XS']))

                  with open(xanesfolder / "xanes_.in", "w") as f:
                      f.write(xinput(mode, iabs[i], dirs,
                                     ortho(atoms.cell, dir - 1,
                                           (dir2 if dir==2 else dir1)-1), xsJSON['XS'],
                                     plot=True))
"""

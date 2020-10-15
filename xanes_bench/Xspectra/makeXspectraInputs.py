# John Vinson, 2020
# based on earlier work by Sencer
"""
Make all the xspectra/qe inputs
"""
from ase.atoms import Atoms
from ase.io import write
import spglib
import pathlib
from os import environ as env
import sys
import numpy as np
from xanes_bench import photonSymm
import json
import xanes_bench.Xspectra
import os

module_path = os.path.dirname(xanes_bench.Xspectra.__file__)
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
            "    xepsilon(1) = %d" % dirs[0],
            "    xepsilon(2) = %d" % dirs[1],
            "    xepsilon(3) = %d" % dirs[2]]

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

    xs_fn = os.path.join(module_path, 'xspectra.json')
    with open (xs_fn, 'r') as fd:
        xsJSON = json.load(fd)

    atoms = smaller( unitCell )

    us = {}
    ph = []
    photonSymm( atoms, us, ph, params['photonOrder'])
    # TODO: Add photonSymm call

    print( ph )

    symm= spglib.get_symmetry((atoms.get_cell(),
                             atoms.get_scaled_positions(),
                             atoms.get_atomic_numbers()),
                             symprec=0.1, angle_tolerance=15)

    equiv = symm['equivalent_atoms']

    symbols = atoms.get_chemical_symbols()


    xsJSON['QE']['electrons']['conv_thr'] = params['defaultConvPerAtom'] * len( symbols )

    sssp_fn = os.path.join(module_path, 'SSSP_precision.json')
    with open (sssp_fn, 'r' ) as pspDatabaseFile:
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


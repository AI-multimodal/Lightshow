# Fanchen Meng, 2022
# John Vinson, 2020
# based on earlier work by Sencer
"""
Make all the xspectra/qe inputs
"""
from ase.atoms import Atoms
from pymatgen.core import Structure
from pymatgen.io.pwscf import PWInput
from pymatgen.io.ase import AseAtomsAdaptor as ase

import xanes_bench.Xspectra
from xanes_bench.photonSym import photonSymm
from xanes_bench.General.kden import printKgrid #, readKgrid, returnKDen, returnKpoint, returnKgridList
from xanes_bench.utils import find_kpts

import json
import numpy as np
from pathlib import Path
import spglib
import os, shutil
import bz2, base64, hashlib
from collections import OrderedDict

module_path = Path(xanes_bench.Xspectra.__path__[0])

def ortho(lat, i, j):
    o = np.linalg.solve(lat.T, np.cross(lat[i], lat[j]))
    return tuple(o / np.linalg.norm(o))

def unpackPsps( ecutwfc, ecutrho, pspDatabaseRoot, DatabaseDir, symbols, folder, needWfn=False ):
    psp = {}
    pspDatabaseName = pspDatabaseRoot + '.json'
    sssp_fn = os.path.join( DatabaseDir, pspDatabaseName)
    with open (sssp_fn, 'r' ) as pspDatabaseFile:
        pspDatabase = json.load( pspDatabaseFile )
    minSymbols = set( symbols )
    for symbol in minSymbols:
        print( symbol )
        print( pspDatabase[ symbol ]['filename'] )
        psp[symbol] = pspDatabase[ symbol ]['filename']
        if ecutwfc < pspDatabase[ symbol ]['cutoff']:
            ecutwfc = pspDatabase[ symbol ]['cutoff']
        if ecutrho < pspDatabase[ symbol ]['rho_cutoff']:
            ecutrho = pspDatabase[ symbol ]['rho_cutoff']
#        if xsJSON['QE']['system']['ecutwfc'] < pspDatabase[ symbol ]['cutoff']:
#            xsJSON['QE']['system']['ecutwfc'] = pspDatabase[ symbol ]['cutoff']
#        if xsJSON['QE']['system']['ecutrho'] < pspDatabase[ symbol ]['rho_cutoff']:
#            xsJSON['QE']['system']['ecutrho'] = pspDatabase[ symbol ]['rho_cutoff']

    pspDatabaseName = pspDatabaseRoot + '_pseudos.json'
    sssp_fn = Path(DatabaseDir)/Path(pspDatabaseName)
    with open (sssp_fn, 'r' ) as p:
        pspJSON = json.load( p )
    for symbol in minSymbols:
        fileName = psp[symbol]
        pspString = bz2.decompress(base64.b64decode( pspJSON[fileName] ))
        print( 'Expected hash:  ' + pspDatabase[symbol]['md5'] )
        print( 'Resultant hash: ' + hashlib.md5( pspString ).hexdigest() )
        with open( folder / ".." / fileName, 'w' ) as f:
            f.write( pspString.decode("utf-8") )

    if needWfn:
        for symbol in minSymbols:
            if 'wfc' not in pspDatabase[ symbol ]:
                print( "WFC not stored corectly in " + pspDatabaseRoot + " for element " + symbol )
                return False
            fileName = pspDatabase[ symbol ]['wfc']
            pspString = bz2.decompress(base64.b64decode( pspJSON[fileName] ))
            print( 'Expected hash:  ' + pspDatabase[symbol]['wfc_md5'] )
            print( 'Resultant hash: ' + hashlib.md5( pspString ).hexdigest() )
            with open( folder / ".." / 'Core.wfc', 'w' ) as f:
                f.write( pspString.decode("utf-8") )

    return psp, ecutwfc, ecutrho


def xinput(mode, iabs, dirs, xkvec, XSparams: dict, plot=False):
    ''' construct input file for XSpectra calculation
        
        Parameters
        ----------
        mode : str, mandatory
            "dipole" or "quadrupole"
        iabs : int, mandatory
            the index of the absorbing element in scf calculation
        dirs : str, mandatory
            TODO
        xkvec : tuple, mandatory
            TODO
        XSparams : dict, mandatory
            paramers parsed to XSpectra calculation
        plot : boolen, optional
            controls xonly_plot

        Returns
        -------
        string of the XSpectra input file
    '''
    inp =  ["&input_xspectra",
            "    calculation = 'xanes_%s'" % mode,
            "    edge = '" + XSparams['input_xspectra']['edge'] + "'",
            "    prefix = 'pwscf'",
            "    outdir = '../'",
            "    xniter = " + str(XSparams['input_xspectra']['xniter']),
            "    xiabs = %d" % iabs,
            "    xerror = " + str(XSparams['input_xspectra']['xerror']),
            "    !wf_collect = .true.",
            "    xcoordcrys = .false.",
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
        # use constant smearing value: 0.89 eV
        # Table IIA in Campbell and Papp (2001) https://doi.org/10.1006/adnd.2000.0848 
        inp += [
            "    gamma_mode = 'constant'",
            "    xgamma = 0.89 ",
            "/"]

    inp += ["&pseudos",
            "    filecore = '../../../Core.wfc'",
            "    r_paw(1) = 1.79", # hard-coded to Ti w/ core-hole
            "/",
            "&cut_occ",
            "    cut_desmooth = " + str( XSparams['cut_occ']['cut_desmooth']),
            "/",
            XSparams['kpts']['kpts'] + " " + XSparams['kpts']['shift'] ]
    return '\n'.join(inp) + '\n'

def makeXspectra( mpid, unitCell: Structure, params: dict ):
    ''' build the input files for the XSpectra workflow.
        
        Parameters
        ----------
        mpid : str, mandatory
            id for structure at Materials Project
        structure : structure in ase.atoms format
        params : dict, mandatory
            TODO
    '''
    # determine element and edge
    element = params['species']
    edge = params['edge']

    # load the default input for XSpectra
    xs_fn = module_path / 'xspectra.json'
    with open (xs_fn, 'r') as fd:
        xsJSON = json.load(fd)
    xsJSON['XS_controls']['element'] = element
    xsJSON['XS_controls']['edge'] = edge
    # target element
    symTarg = xsJSON['XS_controls']['element']
    # build supercell | convert to pymatgen.
    #unitCell = smaller( structure, Rmin=float(xsJSON['XS_controls']['Rmin']) )

    if float(xsJSON['XS_controls']['Rmin']) >= 9:
        # use Gamma point for ground state calculations (es.in and gs.in)
        kpoints = [1,1,1]
    else:
        kpoints = find_kpts(unitCell)

    if params['scf.kpoints'] is not None:
        xs_kpoints = params['scf.kpoints']
    xsJSON['XS']['kpts']['kpts'] = f"{xs_kpoints[0]} {xs_kpoints[1]} {xs_kpoints[2]}"

    us = {}
    symm = spglib.get_symmetry((unitCell.lattice.matrix,
                                unitCell.frac_coords,
                                np.array(unitCell.atomic_numbers)),
                                symprec=0.1, angle_tolerance=15)

    equiv = symm['equivalent_atoms']

    use_photonSymm = True
    ph = []
    if use_photonSymm:
        photonSymm(unitCell, us, ph, params['photonOrder'])
    else:
        directions = {1, 2, 3}
        for dir in range(3):
            dirs = np.zeros(4)
            dirs[-1] = 1.0
            odirs = dirs.copy()
            dirs[dir] = 1.0
            odirs[(dir-1) % 3] = 1.0
            ph.append({
                        "dipole": dirs,
                        "quad":  [odirs]
                      })
        for i in equiv:
            if i in us:
                us[i] = us[i] + 1
            else:
                us[i] = 1
    # get element symbols
    symbols = [str(i).split()[-1] for i in unitCell.species]
    
    xsJSON['QE']['electrons']['conv_thr'] = params['defaultConvPerAtom'] * len( symbols )

    json_dir = params['json_dir']

    folder = Path.cwd() / json_dir / "mp_structures" / mpid / "XS" / \
            Path(f"Spectra-{xs_kpoints[0]}-{xs_kpoints[1]}-{xs_kpoints[2]}")
    folder.mkdir(parents=True, exist_ok=True)

    pspDatabaseRoot = xsJSON['XS_controls']['psp_json']
    DatabaseDir = module_path / '..' / 'pseudos' / 'data' 

    ecutwfc = xsJSON['QE']['system']['ecutwfc']
    ecutrho = xsJSON['QE']['system']['ecutrho']
    psp, ecutwfc, ecutrho = unpackPsps( ecutwfc, ecutrho, pspDatabaseRoot, DatabaseDir, symbols, folder )

    pspDatabaseRoot = xsJSON['XS_controls']['core_psp_json']
    psp2, ecutwfc, ecutrho = unpackPsps( ecutwfc, ecutrho, pspDatabaseRoot, DatabaseDir, 
                                        [xsJSON['XS_controls']['element']], folder, needWfn=True )

    ##TODO for magnetic systems need a more sophisticated system to append numeral
    xsJSON['QE']['system']['ecutwfc'] = ecutwfc
    xsJSON['QE']['system']['ecutrho'] = ecutrho
    xsJSON['QE']['control']['pseudo_dir'] = "../"

    gs_in = PWInput(unitCell, pseudo=psp, control=xsJSON['QE']['control'], 
                    system=xsJSON['QE']['system'], electrons=xsJSON['QE']['electrons'],
                    kpoints_grid=kpoints)
    gs_in.write_file(str(folder / "gs.in"))

    for i in psp2:
        psp[ i + '+' ] = psp2[i]

    psp = OrderedDict(psp)
    for i,j in enumerate(psp.keys()):
        if j == symTarg+'+':
            iabs = i+1
            
    prev = None
    for i, sym in enumerate(symbols):
        if i == equiv[i] and sym == symTarg:
            if prev is not None:
                unitCell[prev] = element

            unitCell[i] = element + '+'
            prev = i

            subfolder = folder / str(i)
            subfolder.mkdir(parents=True, exist_ok=True)
            xsJSON['QE']['control']['pseudo_dir'] = "../../"

            es_in = PWInput(unitCell, pseudo=psp, control=xsJSON['QE']['control'],
                          system=xsJSON['QE']['system'], electrons=xsJSON['QE']['electrons'],
                          kpoints_grid=kpoints)
            es_in.write_file(str(subfolder / "es.in"))
            unitCell[i] = element #'Ti'
            # OCEAN photon labeling is continuous, so we will do that here too
            #  not sure that we will actually want dipole-only spectra(?)
            totalweight = 0
            for photon in ph:
                totalweight += photon["dipole"][3]
          
            photonCount = 0
            for photon in ph:
                photonCount += 1
                dir1 = photon["dipole"][0:3]
                dir2 = dir1
                weight = photon["dipole"][3] * us[i] / totalweight
                mode = "dipole"
                xanesfolder = subfolder / ("%s%d" % (mode, photonCount))
                xanesfolder.mkdir(parents=True, exist_ok=True)
                with open(xanesfolder / "xanes.in", "w") as f:
                    f.write(xinput(mode, iabs, dir1,
                                           dir2, xsJSON['XS']))

                with open(xanesfolder / "xanes_.in", "w") as f:
                    f.write(xinput(mode, iabs, dir1,
                                 dir2,  xsJSON['XS'], plot=True)) 

                with open(xanesfolder / "weight.txt", "w") as f:
                    f.write( str(weight) + "\n" )


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
                        weight = quad[3] * us[i] / totalweight
                        mode = "quadrupole"
                        xanesfolder = subfolder / ("%s%d" % (mode, photonCount))
                        xanesfolder.mkdir(parents=True, exist_ok=True)
                        with open(xanesfolder / "xanes.in", "w") as f:
                            f.write(xinput(mode, iabs[i], dir1,
                                                   dir2, xsJSON['XS']))

                        with open(xanesfolder / "xanes_.in", "w") as f:
                            f.write(xinput(mode, iabs[i], dir1,
                                         dir2,  xsJSON['XS'], plot=True))

                        with open(xanesfolder / "weight.txt", "w") as f:
                            f.write( str(weight) + "\n" )



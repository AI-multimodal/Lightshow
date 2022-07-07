# Fanchen Meng, 2022
# This is a simple tool to find the rmsd E(n,k) between two QE runs
"""
#TODO 
1. move k-point maps (prints and formats) to a separate function or to integers in conjunction
   with the k-point mesh and shift
"""

import xml.etree.ElementTree as ET
import numpy as np
import os
from os import environ as env
from os import path
from scipy.optimize import minimize_scalar

from xanes_bench.EXCITING.parseExciting import readEigval
from xanes_bench.VASP.parseVASP import *


# Hartree to eV based on 2018 CODATA
Ha_c2018 = np.float64( 27.211386245988 )
BroadenParam = 0.01 / Ha_c2018
BroadenParamMetal = 0.3 / Ha_c2018


# Broadening cannot be 0!
def fermiFactor( eig, eFermi, broaden ):
#    ff = ( eig - eFermi )/broaden
    s = (eig - eFermi )/broaden
    if s > 20:
        ff = np.exp( -s ) + np.exp( -2.0*s)
    elif s < -20 :
        ff = 1.0 - np.exp( s ) + np.exp( 2.0*s)
    else:
        ff = 1.0 / ( np.exp( s ) + 1 )
    return ff

# This takes the k-point mesh nk (and eventually its shift kshift)
# And then using the rotation matrixes sym
# returns a complete all-to-all mapping of the k-points kmap
def mapFullKpoints( nk: np.array, kshift: np.array, rots: np.array ):
    xkg = np.zeros((nk[0]*nk[1]*nk[2], 3))
    for i in range(0,int(nk[0])):
        for j in range(0,int(nk[1])):
            for k in range(0,int(nk[2])):
                n = k + j*nk[2] + i*nk[1]*nk[2]
                #JTV check type conversion/rounding
                xkg[n,0] = i/nk[0] # Here we would want to take care of shifted grids
                xkg[n,1] = j/nk[1]
                xkg[n,2] = k/nk[2]

    #  equiv(nk) =nk : k-point nk is not equivalent to any previous k-point
    #  equiv(nk)!=nk : k-point nk is equivalent to k-point equiv(nk)
    equiv = np.zeros((nk[0]*nk[1]*nk[2]),dtype=np.int32)
    for i in range(0,nk[0]*nk[1]*nk[2]):
        equiv[i] = i


    for ik in range(0,nk[0]*nk[1]*nk[2]):
        if equiv[ik] != ik :  # Skip if k-point has already been found to match an earlier one
            continue
        for s in rots:
            xkr = s @ xkg[ik]
            for i in range(0,3):
                xkr[i] = xkr[i] - round( xkr[i] )
            xx = xkr[0]*nk[0]  # Here we would want to take care of shifted grids
            yy = xkr[1]*nk[1]
            zz = xkr[2]*nk[2]
            if abs(xx-round(xx)) <= 0.00001 and abs(yy-round(yy)) <= 0.00001 and abs(zz-round(zz)) <= 0.00001:
                # Here we would want to take care of shifted grids
                i = round( xkr[0]*nk[0] + 2*nk[0] ) % nk[0]
                j = round( xkr[1]*nk[1] + 2*nk[1] ) % nk[1]
                k = round( xkr[2]*nk[2] + 2*nk[2] ) % nk[2]
                n = k + j*nk[2] + i*nk[1]*nk[2]
                if n > ik and equiv[n]==n:
                    equiv[n]=ik
                elif equiv[n]!=ik or n<ik:
                    print( "something wrong in the checking algorithm", n, ik, equiv[n])

            if True:  # time-reversal flag
                xx = -xkr[0]*nk[0]  # Here we would want to take care of shifted grids
                yy = -xkr[1]*nk[1]
                zz = -xkr[2]*nk[2]

                if abs(xx-round(xx)) <= 0.00001 and abs(yy-round(yy)) <= 0.00001 and abs( zz - round(zz) ) <= 0.00001:
                    # Here we would want to take care of shifted grids
                    i = round( -xkr[0]*nk[0] + 2*nk[0] ) % nk[0]
                    j = round( -xkr[1]*nk[1] + 2*nk[1] ) % nk[1]
                    k = round( -xkr[2]*nk[2] + 2*nk[2] ) % nk[2]
                    n = k + j*nk[2] + i*nk[1]*nk[2]
                    if n > ik and equiv[n]==n:
                        equiv[n]=ik
                    elif equiv[n]!=ik or n<ik:
                        print( "something wrong in the checking algorithm")

    totReduced = 0
    for ik in range(0,nk[0]*nk[1]*nk[2]):
        if equiv[ik] == ik:
#            print( xkg[ik] )
            totReduced += 1

#    print( totReduced)

    delta = np.float64( 0.0000000000001 )
    kpointFullMap = {}
    for ik in range(0,nk[0]*nk[1]*nk[2]):
        kvecString = '{:16.8f} {:16.8f} {:16.8f}'.format( xkg[ik][0]+delta, xkg[ik][1]+delta, xkg[ik][2]+delta )
        if equiv[ik] == ik:
            kpointFullMap[kvecString] = []
        else:
            targString = '{:16.8f} {:16.8f} {:16.8f}'.format( xkg[equiv[ik]][0]+delta, xkg[equiv[ik]][1]+delta, xkg[equiv[ik]][2]+delta )
            newList = [targString]
            for i in kpointFullMap[targString]:
                newList.append( i )
                kpointFullMap[ i ].append( kvecString )
            kpointFullMap[targString].append( kvecString )
            kpointFullMap[kvecString] = newList
    return kpointFullMap

def eigPrint( omega, XSkptDict, OkptDict, Okmap ):
    k = 0
    for kpt in XSkptDict:
        k += 1
        okpt = kpt
        if kpt not in OkptDict:
            if kpt not in Okmap:
                print(kpt)
                print(list(Okmap.keys()))
                print( "Incompatible k-point meshes!" )
                exit()
            else:
    # The method used below for the full range is probably better/faster
                for testKpt in OkptDict:
                    for i in Okmap[testKpt]:
    #                    print( i, kpt )
                        if i == kpt:
                            print( "Success: ", testKpt, kpt )
                            okpt = testKpt
                            break
            if okpt == kpt:
                print( "Incompatible k-point meshes!" )
                print( kpt )
                print( Okmap[kpt] )
                print("############")
                print( list(OkptDict.keys()))
                print("############")
                for testKpt in Okmap[kpt]:
                    print( Okmap[testKpt] )
                exit()

        weight = np.float64( XSkptDict[kpt]["weight"] )

        XSeigs = np.asarray( XSkptDict[kpt]["eigenvalues"], dtype=np.float64 )
        Oeigs = np.asarray( OkptDict[okpt]["eigenvalues"], dtype=np.float64 )

        print( "#####" )
        print( kpt )
        print( okpt )
        print( "#####" )
        for j in range( min(len( Oeigs ), len(XSeigs) )  ):
            print( "{:16.8f} {:16.8f} {:16.8f} {:16.8f}".format( XSeigs[j], XSeigs[j]+omega, Oeigs[j], (XSeigs[j]-Oeigs[j]+omega)*Ha_c2018 ) )

        if k > 1:
            break


# Function for calculating the rmsd between two sets of energies
#TODO better error handling?
# omega is the relative offset between them
# XSkptDict and OkptDict are dictionaries containing the two sets of energies

## Optional arguments
# 
def eigRMSD( omega, XSkptDict, OkptDict, Okmap, XSupper=0.0, Oupper=0.0, bUpper=-1.0, 
             XSlower=0.0, Olower=0.0, bLower=-1.0, returnDelta=False ) :
    rmsd = np.float64( 0.0 )
    bandWeight = np.float64( 0.0 )
    maxDelta = np.float64( 0.0 )
    minAllowedWeight = np.float64( 0.00000001 )
    for kpt in XSkptDict:
        okpt = kpt
        if kpt not in OkptDict:
            if kpt not in Okmap:
                print(kpt)
                print(list(Okmap.keys()))
                print( "Incompatible k-point meshes!" )
                exit()
            else:
    # The method used below for the full range is probably better/faster
                for testKpt in OkptDict:
                    for i in Okmap[testKpt]:
    #                    print( i, kpt )
                        if i == kpt:
    #                        print( "Success: ", testKpt, kpt )
                            okpt = testKpt
                            break
            if okpt == kpt:
                print( "Incompatible k-point meshes!" )
                print( kpt )
                print( Okmap[kpt] )
                print("############")
                print( list(OkptDict.keys()))
                print("############")
                for testKpt in Okmap[kpt]:
                    print( Okmap[testKpt] )
                exit()

        weight = np.float64( XSkptDict[kpt]["weight"] )

        XSeigs = np.asarray( XSkptDict[kpt]["eigenvalues"], dtype=np.float64 )
        Oeigs = np.asarray( OkptDict[okpt]["eigenvalues"], dtype=np.float64 )

        for j in range( min(len( Oeigs ), len(XSeigs) )  ):
            if bUpper > 0 :
                off = fermiFactor( Oeigs[j], Oupper, bUpper )
                xsff = fermiFactor( XSeigs[j], XSupper, bUpper )
                ff = np.sqrt( xsff * off )
            else:
                ff = np.float64( 1.0 )
            if bLower > 0 :
                # Note for lower bound just flip the eigenvalue and energy cut-off
                off = fermiFactor( Olower, Oeigs[j], bLower )
                xsff = fermiFactor( XSlower, XSeigs[j], bLower )
                ff = ff * np.sqrt( xsff * off )
            bandWeight += weight * ff
            rmsd += weight*(XSeigs[j]-Oeigs[j]+omega)**2*ff
            m = ff*abs( XSeigs[j]-Oeigs[j]+omega )
            if m > maxDelta:
                maxDelta = m

    #TODO 
    # Need better error handling for this!
    if bandWeight < minAllowedWeight:
        print( "Not enough band weight to calculated rmsd!" )
        rmsd = max( minAllowedWeight*10000, rmsd )
        bandWeight = minAllowedWeight

    if returnDelta:
        return np.sqrt( rmsd/bandWeight ), maxDelta
    return np.sqrt( rmsd/bandWeight )


#def parseQE( filename: str, kmesh: np.array, kshift: np.array, kmap: dict, nelectron: int, energyNK: dict ):
# Opens the QE output xml file filename
# Returns:
# np.arrays for the k-point mesh, k-point shifts
# dictionary of the k-point mapping 
# integer of the number of electrons in the run
# dictionary of the eneries and weights by k-point 
def parseQE( filepath: str ):
    filename = os.path.join( filepath, "pwscf.xml" )
    tree = ET.parse( filename )
    root = tree.getroot()

    # We will want to convert from the QE format to "crystal" format"
    cellList = []
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[0].text.split() )
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[1].text.split() )
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[2].text.split() )
    acell = np.asarray( cellList, dtype=np.double )
    alat = np.float64( root.find('output').find('atomic_structure').get("alat") )

    # Count total electrons, will be useful for figuring you occupied bands
    nelectron = int(float( root.find('output').find('band_structure').find('nelec').text ))
    occBands = int( nelectron / 2 )

    # And grab total number of bands
    nband = int(root.find('output').find('band_structure').find('nbnd').text)

    # Timing 
    time = float(root.find('timing_info').find('total').find('wall').text) 
    print( "N_electron: ", nelectron, "Occupied: ", occBands, "Total bands: ", nband, '{:10.3f} secs'.format(time))

    # rounding fudge to avoid -0
    delta = np.float64( 0.0000000000001 )

    # Grab Fermi
    #TODO deal with insulating case if we make changes
    eFermi = float(root.find('output').find('band_structure').find('fermi_energy').text)
    valMin = np.finfo(np.float64 ).max
    valMax = np.finfo(np.float64 ).min
    condMin = np.finfo(np.float64 ).max
    condMax = np.finfo(np.float64 ).min
    
    energyNK = dict()
    # Now parse and store the first run
    for kpt in root.find('output').find('band_structure').findall('ks_energies'):

        kvec = np.asarray( kpt.find('k_point').text.split(), dtype=np.double )
        kvecReal = acell.dot( kvec ) / alat
        for i in range(0,3):
            if kvecReal[i]+delta < 0.0:
                kvecReal[i]+=1.0
            if kvecReal[i]+delta >= 1.0:
                kvecReal[i]-=1.0
        kvecString = '{:16.8f} {:16.8f} {:16.8f}'.format( kvecReal[0]+delta, kvecReal[1]+delta, kvecReal[2]+delta )
        weight = kpt.find('k_point').get('weight')
      
        eigs = kpt.find('eigenvalues').text.split()

        energyNK[kvecString] = dict()
        energyNK[kvecString]["weight"] = np.float64(weight)/2 #FC
        energyNK[kvecString]["eigenvalues"] = eigs
        for e in np.asarray( eigs, dtype=np.float64 ):
            if e <= eFermi:
                if e > valMax:
                    valMax = e
                if e < valMin:
                    valMin = e
            else:
                if e > condMax:
                    condMax = e
                if e < condMin:
                    condMin = e

    # Grab k-point grid and shift
    kdict = root.find('input').find('k_points_IBZ')[0].attrib
    kmesh = np.array( [kdict['nk1'], kdict['nk2'], kdict['nk3']], dtype=int )
    kshift = np.array( [kdict['k1'], kdict['k2'], kdict['k3']], dtype=int )



    # Now make symmetry matrices 
# FC
#    nrot = int(root.find('output').find('symmetries').find('nrot').text)
#    i = 0
#    sym = np.zeros((nrot,3,3))
#    for s in root.find('output').find('symmetries').findall('symmetry'):
#        sym[i] = np.array( list(map(float,s.find('rotation').text.split())),dtype=int).reshape((3,3),order='F')
#        i += 1
    filename = os.path.join( filepath, "nscf.out" ) 
    with open(filename) as f:
        nscf = f.readlines()
    index = []
    for text in nscf:
        if "Sym. Ops., with inversion, found" in text:
            nrot = int(text.split()[0])
        if "isym" in text:
            index.append(nscf.index(text))
    symm_ = [nscf[ind+2:ind+5] for ind in index]
    sym = []
    for item_lst in symm_:
        tmp = []
        for item in item_lst:
            tmp.append(np.asanyarray(item.split("(  ")[1].split(")")[0].split(), dtype='float'))
        sym.append(np.array(tmp))


    # make kmap dictionary
    kmap = mapFullKpoints( kmesh, kshift, sym )

    bandClips = [ valMin, valMax, condMin, condMax ]

    return kmesh, kshift, kmap, nelectron, energyNK, eFermi, bandClips


def parseEXCITING( filepath: str ):
    # Parse exciting input
    # Get information on kpoints first 
#    kdata_ = np.genfromtxt(os.path.join( "./EXCITING", "KPOINTS.OUT" ),
#            skip_header=1)
    kdata_ = np.genfromtxt(os.path.join( filepath, "KPOINTS.OUT" ),
            skip_header=1)
    kvec_ = np.zeros((kdata_.shape[0],4))
    if kdata_.shape[1] > 6:
        kvec_[:,:] = kdata_[:,[1,2,3,7]] # FC
    else:
        kvec_[:,:] = kdata_[:,[1,2,3,4]] 
    #kvec_[:,:] = kdata_[:,[3,2,1,4]]
    #kvec_[:,:] = kdata_[:,1:5]
    
    ExkptDict = dict()
    for i in range(kvec_.shape[0]):
        kvecString = '{:16.8f} {:16.8f} {:16.8f}'.format( kvec_[i,0],
                kvec_[i,1], kvec_[i,2] )
        ExkptDict[kvecString] = dict()
        ExkptDict[kvecString]["point"] = kvec_[i,:3]
        ExkptDict[kvecString]["weight"] = kvec_[i,3]
    # shift such that the kpoints lie in [-0.5:0.5] interval
    for entry in ExkptDict.keys():
        for i in range(3):
            if ExkptDict[entry]['point'][i] >= 1.0:
                ExkptDict[entry]['point'][i]=ExkptDict[entry]['point'][i]-1.0
            if ExkptDict[entry]['point'][i] <  0.0:
                ExkptDict[entry]['point'][i]=ExkptDict[entry]['point'][i]+1.0
#        print(ExkptDict[entry]['point'])

    # Load up Fermi level
    fname = os.path.join( filepath, "EFERMI.OUT" )
    eFermi = np.genfromtxt( fname ).item()
    valMin = np.finfo(np.float64 ).max
    valMax = np.finfo(np.float64 ).min
    condMin = np.finfo(np.float64 ).max
    condMax = np.finfo(np.float64 ).min
        
    # now read the EIGVAL.OUT file
    fname=os.path.join( filepath, "EIGVAL.OUT")
    eigval_, occBands_exciting, unoccBands_exciting=readEigval(fname)
    print( "Occupied: ", occBands_exciting, "Total bands: ",
            occBands_exciting+unoccBands_exciting)
    for key in ExkptDict:
        ExkptDict[key]["eigenvalues"]=eigval_[key]
        
        for e in np.asarray( eigval_[key], dtype=np.float64 ):
            if e <= eFermi:
                if e > valMax:
                    valMax = e
                if e < valMin:
                    valMin = e
            else:
                if e > condMax:
                    condMax = e
                if e < condMin:
                    condMin = e
        
    bandClips = [ valMin, valMax, condMin, condMax ]

    
    nelectron = occBands_exciting * 2
    return nelectron, ExkptDict, eFermi, bandClips

def parseVASP(filepath):
    fname = os.path.join( filepath, "IBZKPT" )
    kvec, kvec_list, weight = readIBZ(fname)

    ExkptDict = dict()
    for kvecString in kvec.keys():
        ExkptDict[kvecString] = dict()
        ExkptDict[kvecString]["point"] = kvec[kvecString]
        ExkptDict[kvecString]["weight"] = weight[kvecString]
    # shift such that the kpoints lie in [-0.5:0.5] interval
    for entry in ExkptDict.keys():
        for i in range(3):
            if ExkptDict[entry]['point'][i] >= 1.0:
                ExkptDict[entry]['point'][i]=ExkptDict[entry]['point'][i]-1.0
            if ExkptDict[entry]['point'][i] <  0.0:
                ExkptDict[entry]['point'][i]=ExkptDict[entry]['point'][i]+1.0
#        print(ExkptDict[entry]['point'])

    # Load up Fermi level
    fname = os.path.join( filepath, "OUTCAR" )
    eFermi = readFermi(fname)/Ha_c2018
    valMin = np.finfo(np.float64 ).max
    valMax = np.finfo(np.float64 ).min
    condMin = np.finfo(np.float64 ).max
    condMax = np.finfo(np.float64 ).min

    # now read the EIGVAL.OUT file
    fname=os.path.join( filepath, "EIGENVAL")
    eigval_, occBands_vasp, unoccBands_vasp=readEIGEN(fname, kvec_list)
    print( "Occupied: ", occBands_vasp, "Total bands: ",
            occBands_vasp+unoccBands_vasp)
    
    for key in ExkptDict:
        eigval_[key] = eigval_[key]/Ha_c2018
        ExkptDict[key]["eigenvalues"]=eigval_[key]

        for e in np.asarray( eigval_[key], dtype=np.float64 ):
            if e <= eFermi:
                if e > valMax:
                    valMax = e
                if e < valMin:
                    valMin = e
            else:
                if e > condMax:
                    condMax = e
                if e < condMin:
                    condMin = e

    bandClips = [ valMin, valMax, condMin, condMax ]
    print(bandClips)

    nelectron = occBands_vasp * 2
    return nelectron, ExkptDict, eFermi, bandClips

def main():

    mpid = input("Input the mp number: ")
    program = input("[O]cean, [X]spectra, and/or [E]xciting? ")

    AllData = []

    for p in program:
        if( p == 'O' or p == 'o' ):
            fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, 'OCEAN', 
                                    'groundState', "nscf" )
            kmesh, kshift, kmap, nelectron, kptDict, eFermi, clips = parseQE( fileName )
            AllData.append( dict( { "Name" : "OCEAN", "nElectron" : nelectron, 
                                    "kptDict" : kptDict, "eFermi" : eFermi, "clips" : clips, 
                                    "kmesh" : kmesh, "kshift" : kshift, "kmap" : kmap } ) )

        elif( p == 'X' or p == 'x' ):
            fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, 'XS', 
                                    'groundState', "nscf" )
            kmesh, kshift, kmap, nelectron, kptDict, eFermi, clips = parseQE( fileName )
            AllData.append( dict( { "Name" : "XSPECTRA", "nElectron" : nelectron, 
                                    "kptDict" : kptDict, "eFermi" : eFermi, "clips" : clips, 
                                    "kmesh" : kmesh, "kshift" : kshift, "kmap" : kmap } ) )

        elif( p == 'E' or p == 'e' ):
            fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, "EXCITING", "groundState" )
            nelectron, kptDict, eFermi, clips = parseEXCITING( fileName )
            AllData.append( dict( { "Name" : "EXCITING", "nElectron" : nelectron, "kptDict" : kptDict, 
                                    "eFermi" : eFermi, "clips" : clips } ))

        elif( p == 'V' or p == 'v' ):
            fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, "VASP", "groundState" )
            nelectron, kptDict, eFermi, clips = parseVASP( fileName )
            AllData.append( dict( { "Name" : "VASP", "nElectron" : nelectron, "kptDict" : kptDict,
                                    "eFermi" : eFermi, "clips" : clips } ))

        elif( p == 'C' or p == 'c' ):
            customDirName = input("Input the custom directory name: ")
            fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, customDirName,
                                    'groundState', "nscf", "pwscf.xml" )
            if ( path.exists( os.path.join( env['PWD'], "save", "mp_structures", mpid, customDirName,
                                    'groundState', "nscf", "pwscf.xml" ) ) ):
                fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, customDirName,
                                    'groundState', "nscf", "pwscf.xml" )
                kmesh, kshift, kmap, nelectron, kptDict, eFermi, clips = parseQE( fileName )
                AllData.append( dict( { "Name" : customDirName, "nElectron" : nelectron,
                                        "kptDict" : kptDict, "eFermi" : eFermi, "clips" : clips,
                                        "kmesh" : kmesh, "kshift" : kshift, "kmap" : kmap } ) )
            elif ( path.exists( os.path.join( env['PWD'], "save", "mp_structures", mpid, customDirName,
                                    'groundState', "KPOINTS.OUT" ) ) ):
                fileName = os.path.join( env['PWD'], "save", "mp_structures", mpid, customDirName,
                                    'groundState' )
                nelectron, kptDict, eFermi, clips = parseEXCITING( fileName )
                AllData.append( dict( { "Name" : "EXCITING", "nElectron" : nelectron, "kptDict" : kptDict,
                                        "eFermi" : eFermi, "clips" : clips } ) )
            else :
                print("Failed to parse custom options")
        else:
            print("Unrecognized program options")
        

    print( "#          Val min   Val max  Con min  Con max  gap  Con width")
    for i in range( len(AllData) ) :
        print( "{:8s} {: 8.2f} {: 8.2f} {: 8.2f} {: 8.2f} {:7.2f} {:7.2f}"
             .format( AllData[i]["Name"], AllData[i]["clips"][0]*Ha_c2018, AllData[i]["clips"][1]*Ha_c2018,
                      AllData[i]["clips"][2]*Ha_c2018, AllData[i]["clips"][3]*Ha_c2018, 
                      (AllData[i]["clips"][2]-AllData[i]["clips"][1])*Ha_c2018,
                      (AllData[i]["clips"][3]-AllData[i]["clips"][2])*Ha_c2018 ) )



    print( "\nEntire valence band:" )

    for i in range( len(AllData) ) :
        for j in range (i):
    
            print( "  {:8s} vs. {:8s}:".format( AllData[i]["Name"], AllData[j]["Name"] ) )
            if "kmap" in AllData[i]:
                k = i
                m = j
            elif "kmap" in AllData[j]:
                m = i
                k = j
            else :
                print( "!!!At least one data set must have a k-map!!!" )
                exit()

            res = minimize_scalar( eigRMSD, args = (AllData[m]["kptDict"], AllData[k]["kptDict"], 
                    AllData[k]["kmap"], AllData[m]["eFermi"], AllData[k]["eFermi"], BroadenParam), 
                    options={'xtol': 1e-8})
            if res.success:
                print( "    Shift = {: f} eV".format(res.x*Ha_c2018) )
                print( "    RMSD  = {: f} eV".format(res.fun*Ha_c2018) )
            else:
                print( "!!!Optmizing energy shift failed!!!" )
                exit()
            omega = res.x
            rmsd, maxDelta = eigRMSD( omega, AllData[m]["kptDict"], AllData[k]["kptDict"],
                    AllData[k]["kmap"], AllData[m]["eFermi"], AllData[k]["eFermi"], BroadenParam, returnDelta=True)
            print( "    Max D = {: f} eV\n".format(maxDelta*Ha_c2018) )
                


    valenceWindowParam = 20.0

    print( "\nValence band with {:f} eV lower bound:".format(valenceWindowParam))

    for i in range( len(AllData) ) :
        for j in range (i):

            print( "  {:8s} vs. {:8s}:".format( AllData[i]["Name"], AllData[j]["Name"] ) )
            if "kmap" in AllData[i]:
                k = i
                m = j
            elif "kmap" in AllData[j]:
                m = i
                k = j
            else :
                print( "!!!At least one data set must have a k-map!!!" )
                exit()

            lb1 = AllData[m]["clips"][1] - valenceWindowParam/Ha_c2018
            lb2 = AllData[k]["clips"][1]  - valenceWindowParam/Ha_c2018
            res = minimize_scalar( eigRMSD, args = (AllData[m]["kptDict"], AllData[k]["kptDict"],
                    AllData[k]["kmap"], AllData[m]["eFermi"], AllData[k]["eFermi"], BroadenParam, 
                    lb1, lb2, BroadenParam), options={'xtol': 1e-8})

            if res.success:
                print( "    Shift = {: f} eV".format(res.x*Ha_c2018) )
                print( "    RMSD  = {: f} eV".format(res.fun*Ha_c2018) )
            else:
                print( "!!!Optmizing energy shift failed!!!" )
                exit()
            omega = res.x
            rmsd, maxDelta = eigRMSD( omega, AllData[m]["kptDict"], AllData[k]["kptDict"],
                    AllData[k]["kmap"], AllData[m]["eFermi"], AllData[k]["eFermi"], BroadenParam, 
                    lb1, lb2, BroadenParam, True )
            print( "    Max D = {: f} eV\n".format(maxDelta*Ha_c2018) )


    #TODO Might be more efficient to move this up into previous to avoid a second optimization step
    for i in range( len(AllData) ) :
        for j in range (i):

            print( "\n{:8s} vs. {:8s}:".format( AllData[i]["Name"], AllData[j]["Name"] ) )
            if "kmap" in AllData[i]:
                k = i
                m = j
            elif "kmap" in AllData[j]:
                m = i
                k = j
            else :
                print( "!!!At least one data set must have a k-map!!!" )
                exit()

            lb1 = AllData[m]["clips"][1] - valenceWindowParam/Ha_c2018
            lb2 = AllData[k]["clips"][1]  - valenceWindowParam/Ha_c2018
            res = minimize_scalar( eigRMSD, args = (AllData[m]["kptDict"], AllData[k]["kptDict"],
                    AllData[k]["kmap"], AllData[m]["eFermi"], AllData[k]["eFermi"], BroadenParam, 
                    lb1, lb2, BroadenParam), options={'xtol': 1e-8})

            if not res.success:
                print( "!!!Optmizing energy shift failed!!!" )
                exit()
            omega = res.x

            print("  Conduction band with upper bounds\n    Bound  RMSD (eV) Delta (eV)")

            for conductionWindowParam in [ 10.0, 20.0, 30.0, 40.0, 50.0 ]:
                ub1 = AllData[m]["clips"][2] + conductionWindowParam/Ha_c2018
                ub2 = AllData[k]["clips"][2] + conductionWindowParam/Ha_c2018
                rmsd, maxDelta = eigRMSD( omega, AllData[m]["kptDict"], AllData[k]["kptDict"],
                    AllData[k]["kmap"], ub1, ub2, BroadenParam, AllData[m]["eFermi"], AllData[k]["eFermi"],
                    BroadenParam, True )
                print( "    {:5.1f}  {:f}  {:f}".format( conductionWindowParam, rmsd*Ha_c2018, maxDelta*Ha_c2018))


            print("\n  Conduction band with moving 15 eV window")
            print( "    Window RMSD (eV) Delta (eV)")

            for conductionWindowParam in [ 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]:
                XSmid = AllData[m]["clips"][2] + conductionWindowParam/Ha_c2018
                Omid = AllData[k]["clips"][2] + conductionWindowParam/Ha_c2018
                lb1 = max( XSmid - 7.5/Ha_c2018, AllData[m]["eFermi"] )
                lb2 = max( Omid - 7.5/Ha_c2018, AllData[k]["eFermi"] )
                ub1 = XSmid + 7.5/Ha_c2018
                ub2 = Omid + 7.5/Ha_c2018

                rmsd, maxDelta = eigRMSD( omega, AllData[m]["kptDict"], AllData[k]["kptDict"],
                    AllData[k]["kmap"], ub1, ub2, BroadenParam, lb1, lb2, BroadenParam, True )
                print( "    {:5.1f}  {:f}  {:f}".format( conductionWindowParam, rmsd*Ha_c2018, maxDelta*Ha_c2018))

            print( "" )

if __name__ == '__main__':
    main()

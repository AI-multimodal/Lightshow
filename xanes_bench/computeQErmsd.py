# This is a simple tool to find the rmsd E(n,k) between two QE runs
"""
#TODO 
1. move k-point maps (prints and formats) to a separate function or to integers in conjunction
   with the k-point mesh and shift
"""

import xml.etree.ElementTree as ET
import numpy as np
import os
from scipy.optimize import minimize_scalar

from xanes_bench.EXCITING.parseExciting import readEigval


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
def parseQE( filename: str ):
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
    print( "N_electron: ", nelectron, "Occupied: ", occBands, "Total bands: ", nband )

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
        energyNK[kvecString]["weight"] = weight
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
    nrot = int(root.find('output').find('symmetries').find('nrot').text)
    i = 0
    sym = np.zeros((nrot,3,3))
    for s in root.find('output').find('symmetries').findall('symmetry'):
        sym[i] = np.array( list(map(float,s.find('rotation').text.split())),dtype=int).reshape((3,3),order='F')
        i += 1

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
    kvec_[:,:] = kdata_[:,[1,2,3,4]]
#    kvec_[:,:] = kdata_[:,[3,2,1,4]]
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


def main():
    ####
    EXfileName = os.path.join( "./EXCITING", "groundState" )
    EXnelectron, EXkptDict, EXeFermi, EXclips = parseEXCITING( EXfileName )
    print( "EXCITING parsed!")

    #XSfileName = os.path.join( "./XS", "groundState", "pwscf.save", "data-file-schema.xml" )
    XSfileName = os.path.join( "./XS", "groundState", "nscf", "pwscf.xml" )
    XSkmesh, XSkshift, XSkmap, XSnelectron, XSkptDict, XSeFermi, XSclips = parseQE( XSfileName )
    print( "XS parsed.    K-point mesh: {:d} {:d} {:d}.\n"
           .format( XSkmesh[0],XSkmesh[1],XSkmesh[2]))


    #OfileName = os.path.join( "./OCEAN", "groundState", "pwscf.save", "data-file-schema.xml" )
    OfileName = os.path.join( "./OCEAN", "groundState", "nscf", "pwscf.xml" )
    Okmesh, Okshift, Okmap, Onelectron, OkptDict, OeFermi, Oclips = parseQE( OfileName )
    print( "OCEAN parsed. Kpoint mesh {:d} {:d} {:d}.\n"
           .format( Okmesh[0], Okmesh[1], Okmesh[2] ))

    print( "#    Val min   Val max  Con min  Con max  gap  Con width")
    print( "XS {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:7.2f} {:7.2f}".format( XSclips[0]*Ha_c2018, XSclips[1]*Ha_c2018, XSclips[2]*Ha_c2018, XSclips[3]*Ha_c2018, (XSclips[2]-XSclips[1])*Ha_c2018, (XSclips[3]-XSclips[2])*Ha_c2018)  )
    print( " O {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:7.2f} {:7.2f}".format( Oclips[0]*Ha_c2018, Oclips[1]*Ha_c2018, Oclips[2]*Ha_c2018, Oclips[3]*Ha_c2018, (Oclips[2]-Oclips[1])*Ha_c2018, (Oclips[3]-Oclips[2])*Ha_c2018 ) )
    print( "EX {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:7.2f} {:7.2f}".format( EXclips[0]*Ha_c2018, EXclips[1]*Ha_c2018, EXclips[2]*Ha_c2018, EXclips[3]*Ha_c2018, (EXclips[2]-EXclips[1])*Ha_c2018, (EXclips[3]-EXclips[2])*Ha_c2018 ) )


    print( "\nEntire valence band: OCEAN v XS")
    omega = 0
    res = minimize_scalar( eigRMSD, args = (XSkptDict, OkptDict, Okmap, XSeFermi, OeFermi, BroadenParam), options={'xtol': 1e-8})
    if res.success:
        print( "Shift = {:f} eV".format(res.x*Ha_c2018) )
        print( "RMSD  = {:f} eV".format(res.fun*Ha_c2018) )
    else:
        print( "Optmizing energy shift failed" )
        exit()

    omega = res.x
    rmsd, maxDelta = eigRMSD( omega, XSkptDict, OkptDict, Okmap, XSeFermi, OeFermi, BroadenParam, returnDelta=True)
    #print( "RMSD  = {:f} eV".format(rmsd*Ha_c2018) )
    print( "Max D = {:f} eV".format(maxDelta*Ha_c2018) )


    print( "\nEntire valence band: OCEAN v EXCITING")
    omega = 0
    res = minimize_scalar( eigRMSD, args = (EXkptDict, OkptDict, Okmap, EXeFermi, OeFermi, BroadenParam), options={'xtol': 1e-8})
    if res.success:
        print( "Shift = {:f} eV".format(res.x*Ha_c2018) )
        print( "RMSD  = {:f} eV".format(res.fun*Ha_c2018) )
    else:
        print( "Optmizing energy shift failed" )
        exit()

    omega = res.x
    rmsd, maxDelta = eigRMSD( omega, EXkptDict, OkptDict, Okmap, EXeFermi, OeFermi, BroadenParam, returnDelta=True)
    #print( "RMSD  = {:f} eV".format(rmsd*Ha_c2018) )
    print( "Max D = {:f} eV".format(maxDelta*Ha_c2018) )


    valenceWindowParam = 20.0

    print( "\nValence band with {:f} eV lower bound".format(valenceWindowParam))
    omega = 0
    lb1 = XSclips[1] - valenceWindowParam/Ha_c2018
    lb2 = Oclips[1] - valenceWindowParam/Ha_c2018
    res = minimize_scalar( eigRMSD, args = (XSkptDict, OkptDict, Okmap, XSeFermi, OeFermi, BroadenParam,lb1, lb2, BroadenParam), options={'xtol': 1e-8})
    if res.success:
        print( "Shift = {:f} eV".format(res.x*Ha_c2018) )
        print( "RMSD  = {:f} eV".format(res.fun*Ha_c2018) )
    else:
        print( "Optmizing energy shift failed" )
        exit()

    omega = res.x
    rmsd, maxDelta = eigRMSD( omega, XSkptDict, OkptDict, Okmap, XSeFermi, OeFermi, BroadenParam,lb1, lb2, BroadenParam, True)
    #print( "RMSD  = {:f} eV".format(rmsd*Ha_c2018) )
    print( "Max D = {:f} eV".format(maxDelta*Ha_c2018) )



    print( "\nValence band with {:f} eV lower bound OCEAN v EXCITING".format(valenceWindowParam))
    omega2 = 0
    lb1 = EXclips[1] - valenceWindowParam/Ha_c2018
    lb2 = Oclips[1] - valenceWindowParam/Ha_c2018
    res = minimize_scalar( eigRMSD, args = (EXkptDict, OkptDict, Okmap, EXeFermi, OeFermi, BroadenParam,lb1, lb2, BroadenParam), options={'xtol': 1e-8})
    if res.success:
        print( "Shift = {:f} eV".format(res.x*Ha_c2018) )
        print( "RMSD  = {:f} eV".format(res.fun*Ha_c2018) )
    else:
        print( "Optmizing energy shift failed" )
        exit()

    omega2 = res.x
    rmsd, maxDelta = eigRMSD( omega2, EXkptDict, OkptDict, Okmap, EXeFermi, OeFermi, BroadenParam,lb1, lb2, BroadenParam, True)
    #print( "RMSD  = {:f} eV".format(rmsd*Ha_c2018) )
    print( "Max D = {:f} eV".format(maxDelta*Ha_c2018) )

    #eigPrint( omega, EXkptDict, OkptDict, Okmap )
    #exit()


    
    for conductionWindowParam in [ 10.0, 20.0, 30.0, 40.0, 50.0 ]:
        ub1 = XSclips[2] + conductionWindowParam/Ha_c2018
        ub2 = Oclips[2] + conductionWindowParam/Ha_c2018
        rmsd, maxDelta = eigRMSD( omega, XSkptDict, OkptDict, Okmap, ub1, ub2, BroadenParam, 
                                  XSeFermi, OeFermi, BroadenParam, True)
        print( "\nConduction band with {:f} eV upper bound".format(conductionWindowParam))
        print( "RMSD  = {:f} eV".format(rmsd*Ha_c2018) )
        print( "Max D = {:f} eV".format(maxDelta*Ha_c2018) )
        

    print( "\nConduction band with moving 15 eV window")
    print( "Window RMSD (eV) Delta (eV)")
    for conductionWindowParam in [ 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0]:
        XSmid = XSclips[2] + conductionWindowParam/Ha_c2018
        Omid = Oclips[2] + conductionWindowParam/Ha_c2018
        lb1 = max( XSmid - 7.5/Ha_c2018, XSeFermi )
        lb2 = max( Omid - 7.5/Ha_c2018, OeFermi )
        ub1 = XSmid + 7.5/Ha_c2018
        ub2 = Omid + 7.5/Ha_c2018
        rmsd, maxDelta = eigRMSD( omega, XSkptDict, OkptDict, Okmap, ub1, ub2, BroadenParam,
                                  lb1, lb2, BroadenParam, True)
        print( "{:5.1f}  {:f}  {:f}".format(conductionWindowParam, rmsd*Ha_c2018, maxDelta*Ha_c2018))

if __name__ == '__main__':
    main()

# This is a simple tool to find the rmsd E(n,k) between two QE runs
"""
#TODO 
1. move k-point maps (prints and formats) to a separate function or to integers in conjunction
   with the k-point mesh and shift
"""

import xml.etree.ElementTree as ET
import numpy as np
import os

#XSkptDict = dict()
#OkptDict = dict()

# Hartree to eV based on 2018 CODATA
Ha_c2018 = np.float64( 27.211386245988 )


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
    equiv = np.zeros((nk[0]*nk[1]*nk[2]),dtype=np.int)
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
#        print( kvecString )
        

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
    return kmesh, kshift, kmap, nelectron, energyNK
    

XSfileName = os.path.join( "./XS", "groundState", "pwscf.save", "data-file-schema.xml" )
XSkmesh, XSkshift, XSkmap, XSnelectron, XSkptDict = parseQE( XSfileName )
print( "XS parsed.    K-point mesh: {:d} {:d} {:d}. Total bands: {:d}\n"
       .format( XSkmesh[0],XSkmesh[1],XSkmesh[2], XSnelectron))


OfileName = os.path.join( "./OCEAN", "groundState", "pwscf.save", "data-file-schema.xml" )
Okmesh, Okshift, Okmap, Onelectron, OkptDict = parseQE( OfileName )
print( "OCEAN parsed. Kpoint mesh {:d} {:d} {:d}. Total bands: {:d}\n"
       .format( Okmesh[0], Okmesh[1], Okmesh[2], Onelectron))


### Will need to build work on per-code alignment for if XS/Ocean use different semi-core
### and to bring in Exciting
nelectron = XSnelectron
occBands = int( nelectron / 2 )
nband = len( XSkptDict[list(XSkptDict.keys())[0]]["eigenvalues"] )


    
# First calculate the RMSD over the occupied states
rmsd = np.float64( 0.0 )
totalWeight = np.float64( 0.0 )

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
    totalWeight += weight
#    print( weight, totalWeight )

    XSeigs = np.asarray( XSkptDict[kpt]["eigenvalues"], dtype=np.float64 )
    Oeigs = np.asarray( OkptDict[okpt]["eigenvalues"], dtype=np.float64 )

    for j in range( occBands ):
        rmsd += weight*(XSeigs[j]-Oeigs[j])**2

rmsd = np.sqrt( rmsd / totalWeight / occBands )
print( "Occupied :  ", rmsd, rmsd*Ha_c2018)



# Now calculate the RMSD over all the states
rmsd = np.float64( 0.0 )
totalWeight = np.float64( 0.0 )

for kpt in XSkptDict:
    okpt = kpt
    if kpt not in OkptDict:
        if kpt not in Okmap:
            print( "Incompatible k-point meshes!" )
            exit()
        else:
            if kpt not in Okmap:
                print( "Incompatible k-point meshes!" )
                exit()
            for testKpt in Okmap[kpt]:
                if testKpt in OkptDict:
                    okpt = testKpt
                    break
        if okpt == kpt:
            print( "Incompatible k-point meshes!" )
            exit()

    weight = np.float64( XSkptDict[kpt]["weight"] )
    totalWeight += weight

    XSeigs = np.asarray( XSkptDict[kpt]["eigenvalues"], dtype=np.float64 )
    Oeigs = np.asarray( OkptDict[okpt]["eigenvalues"], dtype=np.float64 )

    for j in range( len( Oeigs ) ):
        rmsd += weight*(XSeigs[j]-Oeigs[j])**2

rmsd = np.sqrt( rmsd / totalWeight / nband )
print( "Full range: ", rmsd, rmsd*Ha_c2018 )

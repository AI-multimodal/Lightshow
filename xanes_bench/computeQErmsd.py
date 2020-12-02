# This is a simple tool to find the rmsd E(n,k) between two QE runs

import xml.etree.ElementTree as ET
import numpy as np
import os

from xanes_bench.EXCITING.parseExciting import readEigval

XSkptDict = dict()
OkptDict = dict()
ExkptDict = dict()

# Hartree to eV based on 2018 CODATA
Ha_c2018 = np.float64( 27.211386245988 )


def main():
    tree = ET.parse( os.path.join( "./XS", "pwscf.save", "data-file-schema.xml" ) )
    root = tree.getroot()

    # We will want to convert from the QE format to "crystal" format"
    cellList = []
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[0].text.split() )
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[1].text.split() )
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[2].text.split() )
    acell = np.asarray( cellList, dtype=np.double )
    alat = np.float64( root.find('output').find('atomic_structure').get("alat") )


    # Go ahead and count electrons (occupied bands)
    nelectron = int(float( root.find('output').find('band_structure').find('nelec').text ))
    occBands = int( nelectron / 2 )

    # And grab total number of bands
    nband = int(root.find('output').find('band_structure').find('nbnd').text)
    print( "N_electron: ", nelectron, "Occupied: ", occBands, "Total bands: ", nband )

    delta = np.float64( 0.0000000000001 )

    # Now parse and store the first run
    for kpt in root.find('output').find('band_structure').findall('ks_energies'):
        
        kvec = np.asarray( kpt.find('k_point').text.split(), dtype=np.double )
        kvecReal = acell.dot( kvec ) / alat
        kvecString = '{:18.10f} {:18.10f} {:18.10f}'.format( kvecReal[0]+delta, kvecReal[1]+delta, kvecReal[2]+delta )
        weight = kpt.find('k_point').get('weight')
      
        eigs = kpt.find('eigenvalues').text.split()

        XSkptDict[kvecString] = dict()
        XSkptDict[kvecString]["point"] = kvecReal
        XSkptDict[kvecString]["weight"] = weight
        XSkptDict[kvecString]["eigenvalues"] = eigs
        print( kvecString )


    # Now parse the second run
    tree = ET.parse( os.path.join( "./OCEAN", "pwscf.save", "data-file-schema.xml" ) )
    root = tree.getroot()

    for kpt in root.find('output').find('band_structure').findall('ks_energies'):
        
        kvec = np.asarray( kpt.find('k_point').text.split(), dtype=np.double )
        kvecReal = acell.dot( kvec ) / alat
        kvecString = '{:18.10f} {:18.10f} {:18.10f}'.format( kvecReal[0]+delta, kvecReal[1]+delta, kvecReal[2]+delta )
        weight = kpt.find('k_point').get('weight')
      
        eigs = kpt.find('eigenvalues').text.split()

        OkptDict[kvecString] = dict()
        OkptDict[kvecString]["point"] = kvecReal
        OkptDict[kvecString]["weight"] = weight
        OkptDict[kvecString]["eigenvalues"] = eigs

    # Parse exciting input
    # Get information on kpoints first 
    kdata_ = np.genfromtxt(os.path.join( "./EXCITING", "KPOINTS.OUT" ),
            skip_header=1)
    kvec_ = np.zeros((kdata_.shape[0],4))
    kvec_[:,:] = kdata_[:,[3,2,1,4]]
    #kvec_[:,:] = kdata_[:,1:5]
    
    for i in range(kvec_.shape[0]):
        kvecString = '{:18.10f} {:18.10f} {:18.10f}'.format( kvec_[i,0],
                kvec_[i,1], kvec_[i,2] )
        ExkptDict[kvecString] = dict()
        ExkptDict[kvecString]["point"] = kvec_[i,:3]
        ExkptDict[kvecString]["weight"] = kvec_[i,3]
    # shift such that the kpoints lie in [-0.5:0.5] interval
    for entry in ExkptDict.keys():
        for i in range(3):
            if ExkptDict[entry]['point'][i] > 0.5:
                ExkptDict[entry]['point'][i]=ExkptDict[entry]['point'][i]-1.0
        print(ExkptDict[entry]['point'])
        
    # now read the EIGVAL.OUT file
    fname=os.path.join( "./EXCITING", "EIGVAL.OUT")
    eigval_, occBands_exciting, unoccBands_exciting=readEigval(fname)
    print( "Occupied: ", occBands_exciting, "Total bands: ",
            occBands_exciting+unoccBands_exciting)
    for key in ExkptDict:
        ExkptDict[key]["eigenvalues"]=eigval_[key]
        
    # First calculate the RMSD over the occupied states
    rmsd = np.float64( 0.0 )
    totalWeight = np.float64( 0.0 )

    for kpt in XSkptDict:
        if kpt not in OkptDict:
            print( "Incompatible k-point meshes!" )
            exit()

        weight = np.float64( XSkptDict[kpt]["weight"] )
        totalWeight += weight
    #    print( weight, totalWeight )

        XSeigs = np.asarray( XSkptDict[kpt]["eigenvalues"], dtype=np.float64 )
        Oeigs = np.asarray( OkptDict[kpt]["eigenvalues"], dtype=np.float64 )

        for j in range( occBands ):
            rmsd += weight*(XSeigs[j]-Oeigs[j])**2

    rmsd = np.sqrt( rmsd / totalWeight / occBands )
    print( "Occupied :  ", rmsd, rmsd*Ha_c2018)



    # Now calculate the RMSD over all the states
    rmsd = np.float64( 0.0 )
    totalWeight = np.float64( 0.0 )

    for kpt in XSkptDict:
        if kpt not in OkptDict:
            print( "Incompatible k-point meshes!" )
            exit()

        weight = np.float64( XSkptDict[kpt]["weight"] )
        totalWeight += weight
    #    print( weight, totalWeight )

        XSeigs = np.asarray( XSkptDict[kpt]["eigenvalues"], dtype=np.float64 )
        Oeigs = np.asarray( OkptDict[kpt]["eigenvalues"], dtype=np.float64 )

        for j in range( len( Oeigs ) ):
            rmsd += weight*(XSeigs[j]-Oeigs[j])**2

    rmsd = np.sqrt( rmsd / totalWeight / nband )
    print( "Full range: ", rmsd, rmsd*Ha_c2018 )
if __name__ == '__main__':
    main()

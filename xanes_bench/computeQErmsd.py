# This is a simple tool to find the rmsd E(n,k) between two QE runs

import xml.etree.ElementTree as ET
import numpy as np
import os

XSkptDict = dict()
OkptDict = dict()

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
        OkptDict[kvecString]["weight"] = weight
        OkptDict[kvecString]["eigenvalues"] = eigs


        
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

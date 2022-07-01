# This is a tool to take the data out of QE xml and make data for band structure plots

#TODO 
# 1. Need to figure out labeling high-symmetry points
# 2. Check consistency against EXCITING
# 3. Eventually make part of collect routine?


import xml.etree.ElementTree as ET
import numpy as np
from os.path import exists

# Hartree to eV based on 2018 CODATA
Ha_c2018 = np.float64( 27.211386245988 )


def parseQE( filename: str, outfile ):
    tree = ET.parse( filename )
    root = tree.getroot()

    fermi = np.float64( 0.0 )
    f = root.find('output').find('band_structure').find('fermi_energy')
#    print( f )
    if f is not None:
        fermi = np.float64( f.text )
    print( "fermi = {:f}".format( fermi*Ha_c2018 ) )
    vbm = None

    # We will want to convert from the QE format to "crystal" format"
    cellList = []
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[0].text.split() )
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[1].text.split() )
    cellList.append( root.findall('output')[0].find('atomic_structure').find('cell')[2].text.split() )
    acell = np.asarray( cellList, dtype=np.double )
    alat = np.float64( root.find('output').find('atomic_structure').get("alat") )

    kPathLen = np.float64( 0 )
    allBands = []
    first = True
    distance = []

    for kpt in root.find('output').find('band_structure').findall('ks_energies'):
        kvec = np.asarray( kpt.find('k_point').text.split(), dtype=np.double )
        kvecScaled = kvec * np.pi * 2 / alat

        if first:
            first = False
        else:
            kPathLen += np.linalg.norm( kvecScaled - kvecPrev )

#        print( '{:16.8f} {:16.8f} {:16.8f}'.format( kvec[0], kvec[1], kvec[2] ) )
        eigs = kpt.find('eigenvalues').text.split()
        outfile.write( '{:10.6f}'.format(kPathLen) )
        for e in np.asarray( eigs, dtype=np.float64 ):
            outfile.write( ' {:11.6f}'.format(e*Ha_c2018) )
        outfile.write( "\n" )

        distance.append( kPathLen )
        allBands.append( eigs )
        kvecPrev = kvecScaled

    outfile.close()

    if vbm is None:
        vbm = float(allBands[0][0])
        for i in range(len(allBands)):
            for j in range(len(allBands[i])):
                x = float( allBands[i][j])
                if x <= fermi:
                    if x > vbm:
                        vbm = x

    outfile = open( 'single-column.txt', 'w' )
    for i in range(len(allBands[0])):
        for j in range(len(distance)):
            outfile.write( '{:10.6f} {:11.6f}\n'.format( float(distance[j] ), (float(allBands[j][i])-vbm)*Ha_c2018 ) )
        outfile.write( "\n\n" )
    outfile.close()



    return
    

def main():
    f = "data-file-schema.xml"
    if exists( f ):
        outfile = open( 'test.txt', 'w' )
    else:
        f = "pwscf.xml"
        outfile = open( 'test.txt', 'w' )

    parseQE( f, outfile )


if __name__ == '__main__':
    main()

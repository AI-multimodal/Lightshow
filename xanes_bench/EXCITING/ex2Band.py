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

    # We will want to convert from the QE format to "crystal" format"
    first = True

    distance = []
    allBands = []
    for band in root.findall('band'):
        thisBand = []
        for point in band.findall('point'):
            if first:
                distance.append( point.attrib['distance'] )
            thisBand.append( point.attrib['eval'] )
        first = False
        allBands.append( thisBand )


    tpBands = np.transpose( np.asarray( allBands, dtype=np.float64 ) )

    for j in range(len(distance)):
        outfile.write( '{:10.6f}'.format( float(distance[j] )) )
        for e in tpBands[j]:
            outfile.write( ' {:11.6f}'.format(e*Ha_c2018) )
        outfile.write( "\n" )

    outfile.close()

    vbm = tpBands[0][0]
    for i in range(len(tpBands)):
        for j in range(len(tpBands[i])):
            if tpBands[i][j] <= 0.0:
                if tpBands[i][j] > vbm:
                    vbm = tpBands[i][j]

    outfile = open( 'single-column.txt', 'w' )
    for i in range(len(tpBands[0])):
        for j in range(len(distance)):
            outfile.write( '{:10.6f} {:11.6f}\n'.format( float(distance[j] ), (tpBands[j][i]-vbm)*Ha_c2018 ) )
        outfile.write( "\n\n" )
    outfile.close()
        
            

    return
    

def main():
    f = "bandstructure.xml"
    if exists( f ):
        outfile = open( 'test.txt', 'w' )
    else:
        print( "bandstructure.xml not found!" )
        exit()

    parseQE( f, outfile )


if __name__ == '__main__':
    main()

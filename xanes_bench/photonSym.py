# John Vinson, 2020
"""
Determines set of photon polarization/q-direction along with symmetry unique atoms
"""
#from ase.atoms import Atoms
from pymatgen.core import Structure
import spglib
import pprint
import numpy as np


def BuildOh( n, oh, norm ):
    w = 1
    norm = 6
    if n > 6:
        w = 40
    oh.append( [1,0,0,w] )
    oh.append( [0,1,0,w] )
    oh.append( [0,0,1,w] )
    oh.append( [-1,0,0,w] )
    oh.append( [0,-1,0,w] )
    oh.append( [0,0,-1,w] )

    if n < 6:
        return
    norm = 30

    w = 2
    if n > 6:
        w = 32
    oh.append( [1,1,0,w] )
    oh.append( [1,-1,0,w] )
    oh.append( [-1,1,0,w] )
    oh.append( [-1,-1,0,w] )
    oh.append( [1,0,1,w] )
    oh.append( [1,0,-1,w] )
    oh.append( [-1,0,1,w] )
    oh.append( [-1,0,-1,w] )
    oh.append( [0,1,1,w] )
    oh.append( [0,1,-1,w] )
    oh.append( [0,-1,1,w] )
    oh.append( [0,-1,-1,w] )

    if n < 8:
      return
    norm = 840

    w = 27
    oh.append( [1,1,1,w] )
    oh.append( [-1,1,1,w] )
    oh.append( [1,-1,1,w] )
    oh.append( [1,1,-1,w] )
    oh.append( [-1,-1,1,w] )
    oh.append( [1,-1,-1,w] )
    oh.append( [-1,1,-1,w] )
    oh.append( [-1,-1,-1,w] )





#def photonSym( mode: str, atoms: Atoms, uniqueSites: dict, photons: list ):
def photonSymm( atoms: Structure, uniqueSites: dict, photons: list, order=4 ):
    # for now mode is assumed to be quad

#    uniqueSites = {}
#    photons = []

    symm= spglib.get_symmetry((atoms.lattice.matrix,
                               atoms.frac_coords,
                               np.array(atoms.atomic_numbers)),
                               symprec=0.1, angle_tolerance=15)

    equiv = symm['equivalent_atoms']

    for i in equiv:
        if i in uniqueSites:
            uniqueSites[i] = uniqueSites[i] + 1
        else:
            uniqueSites[i] = 1


    print( uniqueSites )
    photons.append( { "dipole": [1,0,0,1] } )
    photons.append( { "dipole": [0,1,0,1] } )
    photons.append( { "dipole": [0,0,1,1] } )
    return


    ohList = []
    norm = 0
    BuildOh( order, ohList, norm )
    oh1 = np.vstack( ohList )


    finalDipoleList = []
    weight = 0

    symm['rotations'] = np.reshape( np.append(symm['rotations'], [[-1,0,0],[0,-1,0],[0,0,-1]] ), (-1,3,3) )

#    for rot in symm['rotations']:
#        print( rot)


    for sym in oh1:
#        print( sym )
#        trimSym = np.ndarray((3,), buffer=sym, dtype='intc' )
        trimSym = sym[0:3]
        found = False
#        print( "=========------------" )
#        print( trimSym )
#        exit()
        for i, test in enumerate( finalDipoleList ):
#            testSym = np.ndarray((3,), buffer=test, dtype='intc' )
            testSym = test[0:3]
#            print( "------------" )
#            print( testSym, trimSym )
            for rot in symm['rotations']:
#                print( rot)
#                rotSym = rot * trimSym 
                rotSym = rot.dot( trimSym )
#                print( rotSym )
                if( np.all(rotSym == testSym )):
                    found = True
                    weight = sym[3]
                    break
            if found:
                break

        if not found:
            finalDipoleList.append( sym )
        else:
            finalDipoleList[i][3] += weight


    print( "Final dipoles" )
    for test in finalDipoleList:
        print( test )



    for i, dip in enumerate( finalDipoleList ):
        dipVec = dip[0:3]
        weight = dip[3]

        test = np.array( [1,0,0] )
        if( np.all( test == dipVec ) ):
            test = np.array( [0,1,0] )

        quad1 = np.cross( dipVec, test )
        quad2 = np.cross( dipVec, quad1 )
        quad = np.append( quad1, quad2 )
        quad = np.append( quad, np.add( quad1, quad2 ) )
        quad = np.reshape( np.append( quad, np.subtract( quad1, quad2 ) ), (4,3))

        tempQuadList = []
        for q in quad:
            found = False
            for i, test in enumerate( tempQuadList ):
                testSym = test[0:3]
                for rot in symm['rotations']:
                    rotSym = rot.dot( trimSym )
                    if( np.all(rotSym == testSym )):
                        found = True
                        break
                if found:
                    break
            if not found:
                storeQuad = np.append( q, weight )
                tempQuadList.append( storeQuad )
            else:
                tempQuadList[i][3] += weight

#        for q in tempQuadList:
#            print(q)
        photons.append( { "dipole": dip, "quad": tempQuadList, "norm": norm } )

        print( { "dipole": dip, "quad": tempQuadList } )
#    print( photons )

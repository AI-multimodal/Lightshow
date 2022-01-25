""" Prints list of k-point grids 
"""

#from ase.atoms import Atoms
from pymatgen.core import Structure
from ase.units import Bohr
from math import pi
import numpy as np

def printKgrid( unitC: Structure, folder: str ):

    fd = open (str(folder / 'k.txt'), 'w') 
#    recip = atoms.cell.reciprocal().lengths()

    #vol = np.dot( np.cross(  unitC.cell[0],  unitC.cell[1] ), unitC.cell[2] )
    
    #recip = [ np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[1] )/vol ), 
    #          np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[0] )/vol ), 
    #          np.linalg.norm( np.cross(  unitC.cell[0],  unitC.cell[1] )/vol ) ]
    recip = [np.linalg.norm(unitC.lattice.reciprocal_lattice.matrix[i]/2/pi) for i in range(3)]

    klen = 1.0 / recip[0]
    for i in range(3):
        t = 1.0 / recip[i]
        if t < klen:
            klen = t

    kpt = [1,1,1]
#    print( " kx   ky   kz          kden   (inv Bohr)      (inv A)  radius (Bohr)" )
#    print( "{:3d}  {:3d}  {:3d}  {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format( kpt[0], kpt[1], kpt[2], klen/(2*pi), 2*pi*Bohr/klen, 2*pi/klen, klen/Bohr ) )
    fd.write( "#kx   ky   kz     (inv Bohr)     (inv A)   radius (Bohr)  radius(A)\n" )
    fd.write( "{:3d}  {:3d}  {:3d}  {:12.5f} {:12.5f} {:12.5f} {:12.5f}\n".format( kpt[0], kpt[1], kpt[2], 2*pi*Bohr/klen, 2*pi/klen, klen/Bohr, klen ) )

    while klen < 53:
        klen *= 10
        t = kpt[0]/recip[0] * 1.01
        if( t > kpt[1]/recip[1] * 1.01 ):
            t = kpt[1]/recip[1] * 1.01 
        if( t > kpt[2]/recip[2] * 1.01 ):
            t = kpt[2]/recip[2] * 1.01

        for i in range(3):
            kpt[i] = int( t * recip[i] ) + 1
            if klen > kpt[i]/recip[i] :
                klen = kpt[i]/recip[i]
#        print( "{:3d}  {:3d}  {:3d}  {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format( kpt[0], kpt[1], kpt[2], klen/(2*pi), 2*pi*Bohr/klen, 2*pi/klen, klen/Bohr ) )
#        fd.write( "{:3d}  {:3d}  {:3d}  {:12.5f} {:12.5f} {:12.5f} {:12.5f}\n".format( kpt[0], kpt[1], kpt[2], klen/(2*pi), 2*pi*Bohr/klen, 2*pi/klen, klen/Bohr ) )
        fd.write( "{:3d}  {:3d}  {:3d}  {:12.5f} {:12.5f} {:12.5f} {:12.5f}\n".format( kpt[0], kpt[1], kpt[2], 2*pi*Bohr/klen, 2*pi/klen, klen/Bohr, klen ) )


    fd.close()


def returnKgridList( unitC: Structure, maxLen=53.0 ):
    kglist = []

    #vol = np.dot( np.cross(  unitC.cell[0],  unitC.cell[1] ), unitC.cell[2] )
    #recip = [ np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[1] )/vol ),
    #          np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[0] )/vol ),
    #          np.linalg.norm( np.cross(  unitC.cell[0],  unitC.cell[1] )/vol ) ]
    recip = [np.linalg.norm(unitC.lattice.reciprocal_lattice.matrix[i]/2/pi) for i in range(3)]

    klen = 1.0 / recip[0]
    for i in range(3):
        t = 1.0 / recip[i]
        if t < klen:
            klen = t

    kpt = [1,1,1,klen]
    kglist.append( kpt.copy() )

    while klen < maxLen:
        klen *= 10
        t = kpt[0]/recip[0] * 1.01
        if( t > kpt[1]/recip[1] * 1.01 ):
            t = kpt[1]/recip[1] * 1.01
        if( t > kpt[2]/recip[2] * 1.01 ):
            t = kpt[2]/recip[2] * 1.01

        for i in range(3):
            kpt[i] = int( t * recip[i] ) + 1
            if klen > kpt[i]/recip[i] :
                klen = kpt[i]/recip[i]
        for i in range(3):
            kpt[i] = int( t * recip[i] ) + 1
            if klen > kpt[i]/recip[i] :
                klen = kpt[i]/recip[i]

        kpt[3] = klen
        kglist.append( kpt.copy() )

    for i in kglist:
        print(i)

    return kglist

def readKgrid( folder: str ):
    klist={}
    with open (str(folder / 'k.txt'), 'r') as fd:
        for line in fd.readlines():
            kx,ky,kz= line.split()[0:3]
            if not kx.isnumeric():
                continue
            else:
                klist[(int(kx),int(ky),int(kz))]=float(line.split()[-1])
    return klist

def returnKDen( unitC: Structure, kpoint ):

    #vol = np.dot( np.cross(  unitC.cell[0],  unitC.cell[1] ), unitC.cell[2] )

    #recip = [ np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[1] )/vol ),
    #          np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[0] )/vol ),
    #          np.linalg.norm( np.cross(  unitC.cell[0],  unitC.cell[1] )/vol ) ]
    recip = [np.linalg.norm(unitC.lattice.reciprocal_lattice.matrix[i]/2/pi) for i in range(3)]

    klen = kpoint[0] / recip[0]
    for i in range(3):
        t = kpoint[i] / recip[i]
        if t < klen:
            klen = t

    return klen


def returnKpoint( unitC: Structure, klen ):
    #vol = np.dot( np.cross(  unitC.cell[0],  unitC.cell[1] ), unitC.cell[2] )

    #recip = [ np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[1] )/vol ),
    #          np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[0] )/vol ),
    #          np.linalg.norm( np.cross(  unitC.cell[0],  unitC.cell[1] )/vol ) ]
    recip = [np.linalg.norm(unitC.lattice.reciprocal_lattice.matrix[i]/2/pi) for i in range(3)]
    kpoint = []
    for i in range(3):
        kpoint.append( int( klen * recip[i] ) + 1 )

    return kpoint

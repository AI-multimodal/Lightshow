""" Prints list of k-point grids 
"""

from ase.atoms import Atoms
from ase.units import Bohr
from math import pi
import numpy as np

def printKgrid( unitC: Atoms, folder: str ):

    fd = open (str(folder / 'k.txt'), 'w') 
#    recip = atoms.cell.reciprocal().lengths()

    vol = np.dot( np.cross(  unitC.cell[0],  unitC.cell[1] ), unitC.cell[2] )
    
    recip = [ np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[1] )/vol ), 
              np.linalg.norm( np.cross(  unitC.cell[2],  unitC.cell[0] )/vol ), 
              np.linalg.norm( np.cross(  unitC.cell[0],  unitC.cell[1] )/vol ) ]

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

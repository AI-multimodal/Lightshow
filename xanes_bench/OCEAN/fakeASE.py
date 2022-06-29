# coding: latin-1
# J Vinson, 2020

""" Reads/write OCEAN files
"""
#from ase.atoms import Atoms
from pymatgen.core import Structure
#from collections import OrderedDict
import os
import sys
#from pymatgen.io.vasp.sets import MPStaticSet
from ase.units import Bohr
import numpy as np

def write_ocean_in( filename: str, structure: Structure, input_data: dict ):

    filename = os.path.expanduser(filename)
    mode = 'w'
    fd = open(filename, mode)

    #Change screen.nkpt -Int to triplet of ints
    if 'screen.nkpt' in input_data and len(input_data['screen.nkpt'].split()) == 1:
       input_data['screen.nkpt'] = oceanKptSampling( structure.lattice.matrix, float(input_data['screen.nkpt']) )


    input_data_str = []
    for key in input_data:
      input_data_str.append( str(key) + ' { ' + str(input_data[key]) + ' }\n' )

    fd.write( ''.join(input_data_str ))
    # TODO: corresponding pymatgen methods? using StaticSet?
    #if any(atoms.get_initial_magnetic_moments()):
    #    raise NameError( 'Spin=2 not implemented yet' )
    # TODO: how to choose the spin for OCEAN?
#    static = MPStaticSet(structure)
#    if any(static.incar['MAGMOM']) > 0 or any(static.incar['MAGMOM']) < 0 :
#        raise NameError( 'Spin=2 not implemented yet' )

#    atomic_species = OrderedDict()
#    atomic_species_str = []    
#    znucl = []
#    typat = []
    
#    ispecies = 0

#    for atom in atoms:
#        if atom.symbol not in atomic_species:
#            ispecies += 1
#            atomic_species[atom.symbol] = ispecies  # just a placeholder
#  
#        znucl.append( str( atom.numbers ) + ' ' )
#        typat.append( str( atomic_species[atom.symbol] ) + ' ' )
#        atomic_positions_str.append( '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}\n'.format(
#              coords=[atom.a, atom.b, atom.c] ))


    #species = sorted(set(atoms.numbers))
    species = sorted(set(structure.atomic_numbers))

    fd.write('znucl {{ {} }}\n'.format(' '.join(str(Z) for Z in species)))
    fd.write('typat')
    fd.write('{\n')
    types = []
    for Z in structure.atomic_numbers:
        for n, Zs in enumerate(species):
            if Z == Zs:
                types.append(n + 1)
    n_entries_int = 20  # integer entries per line
    for n, type in enumerate(types):
        fd.write(' %d' % (type))
        if n > 1 and ((n % n_entries_int) == 1):
            fd.write('\n')
    fd.write(' }\n')

    atomic_positions_str = []
    for atom in structure:
        atomic_positions_str.append( '{coords[0]:.10f} {coords[1]:.10f} {coords[2]:.10f}\n'.format(
              coords=[atom.a, atom.b, atom.c] ))

    fd.write( 'xred {\n' )
    fd.write( ''.join(atomic_positions_str))
    fd.write( '}\n' )

    fd.write( 'acell {{ {acell[0]} {acell[0]} {acell[0]} }} \n'.format( acell=[1/Bohr] ) )
    
    fd.write( 'rprim {{ {cell[0][0]:.14f} {cell[0][1]:.14f} {cell[0][2]:.14f}\n'
              '        {cell[1][0]:.14f} {cell[1][1]:.14f} {cell[1][2]:.14f}\n'
              '        {cell[2][0]:.14f} {cell[2][1]:.14f} {cell[2][2]:.14f}  }}\n'
                   ''.format(cell=structure.lattice.matrix))



def oceanKptSampling( cell, kpt ):
    # Should normalize each cell dim by Bohr, or we can just do it the once for volume
    v = abs( np.dot(np.cross(cell[0], cell[1]), cell[2] ))/Bohr
    b1 = 2*np.pi*np.linalg.norm(np.cross(cell[1], cell[2]))/v
    b2 = 2*np.pi*np.linalg.norm(np.cross(cell[0], cell[2]))/v
    b3 = 2*np.pi*np.linalg.norm(np.cross(cell[0], cell[1]))/v

    k1 = int( -kpt * b1 ) + 1
    k2 = int( -kpt * b2 ) + 1
    k3 = int( -kpt * b3 ) + 1

    return "{:d} {:d} {:d}".format( k1, k2, k3 )

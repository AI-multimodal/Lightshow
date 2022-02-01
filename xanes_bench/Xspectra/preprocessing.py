from pathlib import Path
from ase.atoms import Atoms
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor as ase
import xanes_bench
import json
import numpy as np
import spglib
def smaller(atoms: Atoms, Rmin=9.0):
    ''' Build the supercell

        Parameters
        ----------
        atoms : structure in ase.atom format
        Rmin : float, optional
            the desired minimum lattice constant

        Returns
        -------
        pymatgen.core.Structure
    '''
    prim =  atoms * ((Rmin / np.linalg.norm(atoms.cell, axis=1)).astype(int) + 1)
    lat, pos, Z = spglib.standardize_cell((atoms.get_cell(),
                                           atoms.get_scaled_positions(),
                                           atoms.get_atomic_numbers()))
    conv = Atoms(Z, cell=lat, positions=pos@lat, pbc=True)
    conv = conv * ((Rmin / np.linalg.norm(conv.cell, axis=1)).astype(int) + 1)
    return ase.get_structure(conv) if len(conv) <= len(prim) else ase.get_structure(prim)

def build_supercell(atoms: Atoms, Rmin=9.0):
    module_path = Path(xanes_bench.Xspectra.__path__[0])
    xs_fn = module_path / 'xspectra.json'
    with open (xs_fn, 'r') as fd:
        xsJSON = json.load(fd)
    structure = smaller( atoms, Rmin=float(xsJSON['XS_controls']['Rmin']) )
    return structure

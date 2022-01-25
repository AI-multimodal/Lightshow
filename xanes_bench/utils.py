import numpy as np
from math import pi
# Return a guess at the number of conduction bands that a given unit-cell volume needs to 
# cover a given energy range (in Ryd)
#JTV
#TODO needs to be unified with groundState.py, both call the same function 
#     and the prefactor should always be the same
def getCondBands( volume, eRange):
    return round( 0.256 * volume * ( eRange**(3/2) ) )

def dielectricGuess( gap ):
    if gap < 0.00001:
        return 1000000
    return 1.0/(2*atan(.164475*gap)/pi)

def find_kpts(cell, cutoff=24.0, max_radii=50.0):
    """Customizes the KPOINTS file.
    Parameters
    ----------
    supercell : paymatgen.core.structure.Structure
        Previously generated supercell.
    cutoff : float, optional
        Description
    max_radii : float, optional
        Description
    Returns
    -------
    pymatgen.io.vasp.Kpoints
    """

    klist = dict()
    rlatt = np.array(cell.lattice.reciprocal_lattice.abc)

    for xx in np.arange(0, 10, 0.2):
        div = np.floor(xx * rlatt) + 1
        divlatt = 2.0 * pi / rlatt * div

        radi = min(divlatt)

        if radi > max_radii:
            break
        else:
            div = tuple(div.astype(int))
            if div not in klist:
                klist[div] = radi

    for key, value in klist.items():
        if value > cutoff:
            k = key
            break

    return  k 


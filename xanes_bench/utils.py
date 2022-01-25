# Fanchen Meng 2022
# based on previous work by J Vinson
import numpy as np
from math import pi
import xanes_bench

from pymatgen.ext.matproj import MPRester

def get_structure(mpid):
    ''' get material strcutre from Materials Project
        
        Parameters
        ----------
        mpid : str, mandatory
            material id in Materials Project

        Returns
        -------
        st : pymatgen.core.Structure
            structure data
        st_dict : dict
            metadata
    '''
    mpkey_fn = Path(xanes_bench.__path__[0]) / "mp.key"
    with open(mpkey_fn, 'r' ) as f:
        mpkey = f.read()
        mpkey = mpkey.strip()

    mpr = MPRester( str(mpkey) )

    try:
        st = mpr.get_structure_by_material_id(mpid, conventional_unit_cell=False)
    except Exception as e:
        print(e)
        print( "Failed to 'get_structure_by_material_id'\nStopping\n" )
        exit()

    st_dict = st.as_dict().copy()
    st_dict["download_at"] = time.ctime()
    try:
        st_dict["created_at"] = mp.get_doc(mpid)["created_at"]
    except Exception as e:
        print(e)
        print( "Failed to 'get_doc'\nStopping\n")
        exit()
    return st, st_dict

# Return a guess at the number of conduction bands that a given unit-cell volume needs to 
# cover a given energy range (in Ryd)
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


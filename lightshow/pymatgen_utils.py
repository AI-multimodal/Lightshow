import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor as ase
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def make_supercell(prim, cutoff=9.0):
    """Used to generate supercell with desired lattice vector: default cutoff
    9 Angstrom, which is am empirical value according to our test. It will
    choose either to use primitive unit cell or to use  conventional unit cell
    depending on which has less atoms.

    Parameters
    ----------
    prim : pymatgen.core.structure.Structure
        The primitive structure.
    cutoff : float, optional
        The supercell cutoff.

    Returns
    -------
    pymatgen.core.structure.Structure
    """

    # If no MAGMOM data, setup for spin polarized calculations
    if "magmom" not in prim.site_properties:
        magmom = [1 for i in range(len(prim))]
        prim.add_site_property("magmom", magmom)
    prim_magmom = prim.site_properties["magmom"]

    latt_prim = np.array(prim.lattice.abc)

    shape_prim = np.ceil(cutoff / latt_prim).astype(int)

    prim_super = ase.get_atoms(prim, magmoms=prim_magmom) * shape_prim

    magmom = prim_super.get_initial_magnetic_moments()

    structure = ase.get_structure(prim_super)

    for idx, site in enumerate(structure):
        site.properties["magmom"] = magmom[idx]

    return structure


def get_inequivalent_site_info(structure):
    """Gets the symmetrically inequivalent sites as found by the
    SpacegroupAnalyzer class from Pymatgen.

    Parameters
    ----------
    structure : pymatgen.core.structure.Structure
        The Pymatgen structure of interest.

    Returns
    -------
    dict
        A dictionary containing three lists, one of the inequivalent sites, one
        for the atom types they correspond to and the last for the multiplicity.
    """

    # Get the symmetrically inequivalent indexes
    inequivalent_sites = (
        SpacegroupAnalyzer(structure)
        .get_symmetrized_structure()
        .equivalent_indices
    )

    # Equivalent indexes must all share the same atom type
    multiplicities = [len(xx) for xx in inequivalent_sites]
    inequivalent_sites = [xx[0] for xx in inequivalent_sites]
    species = [str(structure[xx].specie) for xx in inequivalent_sites]

    return {
        "sites": inequivalent_sites,
        "species": species,
        "multiplicities": multiplicities,
    }

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
    species = [structure[xx].specie.symbol for xx in inequivalent_sites]

    return {
        "sites": inequivalent_sites,
        "species": species,
        "multiplicities": multiplicities,
    }


def get_supercell_indexes_matching_primitive(prim, sc, compare, r):
    info_prim = get_inequivalent_site_info(prim)
    info_sc = get_inequivalent_site_info(sc)

    # These should be the "master indexes" that we're matching against
    primitive_indexes = info_prim["sites"]
    primitive_species = info_prim["species"]

    # These are the supercell inequivalent indexes
    sc_indexes = info_sc["sites"]
    sc_species = info_sc["species"]

    # Get the neighbors of the primitive and supercell inequivalent indexes
    primitive_neighbors = [
        sorted([jj.nn_distance for jj in prim.get_neighbors(prim[ii], r=r)])
        for ii in primitive_indexes
    ]
    # primitive_species = [prim[ii].specie.symbol for ii in primitive_indexes]
    sc_neighbors = [
        sorted([jj.nn_distance for jj in sc.get_neighbors(sc[ii], r=r)])
        for ii in sc_indexes
    ]
    # sc_species = [sc[ii].specie.symbol for ii in sc_indexes]

    # Match
    site_matching = {}
    for cc, (primitive_index, p_species) in enumerate(
        zip(primitive_indexes, primitive_species)
    ):
        primitive_neighbor_distances = primitive_neighbors[cc]
        L_prim = len(primitive_neighbor_distances)

        sc_neighbor_vector_lengths = [len(xx) for xx in sc_neighbors]
        sc_indexes_with_same_length_as_prim = [
            ii
            for ii, (LL, symbol) in enumerate(
                zip(sc_neighbor_vector_lengths, sc_species)
            )
            if LL == L_prim and p_species == symbol
        ]

        # Out of these candidates, we compare the distances
        sc_distances_matches = [
            np.allclose(
                primitive_neighbor_distances[:compare],
                sc_neighbors[ii][:compare],
            )
            for ii in sc_indexes_with_same_length_as_prim
        ]
        # print(sc_indexes_with_same_length_as_prim)
        # print(sc_distances_matches)
        index = sc_distances_matches.index(True)
        # print(sc_indexes_with_same_length_as_prim[index])
        # print(index)

        site_matching[primitive_index] = sc_indexes[
            sc_indexes_with_same_length_as_prim[index]
        ]
    # print(primitive_species)
    return site_matching


def atom_in_structure(atom_symbol, structure):
    """Checks the provided structure to see if the atom_symbol is present
    in it.

    Parameters
    ----------
    atom_symbol : str
    structure : pymatgen.core.structure.Structure

    Returns
    -------
    bool
    """

    return atom_symbol in [s.specie.symbol for s in structure]

from pathlib import Path
from functools import lru_cache

import numpy as np
from pymatgen.core.structure import IStructure
from pymatgen.io.exciting import ExcitingInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# 6 is the default precision used in FEFF by PMG!!!
# See PR I'm going to open soon
ROUNDING_PRECISON = 6

ATOL = 1e-4


ATOMIC_NUMBERS = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "Ne": 10,
    "Na": 11,
    "Mg": 12,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Sc": 21,
    "Ti": 22,
    "V": 23,
    "Cr": 24,
    "Mn": 25,
    "Fe": 26,
    "Co": 27,
    "Ni": 28,
    "Cu": 29,
    "Zn": 30,
    "Ga": 31,
    "Ge": 32,
    "As": 33,
    "Se": 34,
    "Br": 35,
    "Kr": 36,
    "Rb": 37,
    "Sr": 38,
    "Y": 39,
    "Zr": 40,
    "Nb": 41,
    "Mo": 42,
    "Tc": 43,
    "Ru": 44,
    "Rh": 45,
    "Pd": 46,
    "Ag": 47,
    "Cd": 48,
    "In": 49,
    "Sn": 50,
    "Sb": 51,
    "Te": 52,
    "I": 53,
    "Xe": 54,
    "Cs": 55,
    "Ba": 56,
    "La": 57,
    "Ce": 58,
    "Pr": 59,
    "Nd": 60,
    "Pm": 61,
    "Sm": 62,
    "Eu": 63,
    "Gd": 64,
    "Tb": 65,
    "Dy": 66,
    "Ho": 67,
    "Er": 68,
    "Tm": 69,
    "Yb": 70,
    "Lu": 71,
    "Hf": 72,
    "Ta": 73,
    "W": 74,
    "Re": 75,
    "Os": 76,
    "Ir": 77,
    "Pt": 78,
    "Au": 79,
    "Hg": 80,
    "Tl": 81,
    "Pb": 82,
    "Bi": 83,
    "Po": 84,
    "At": 85,
    "Rn": 86,
    "Fr": 87,
    "Ra": 88,
    "Ac": 89,
    "Th": 90,
    "Pa": 91,
    "U": 92,
    "Np": 93,
    "Pu": 94,
    "Am": 95,
    "Cm": 96,
    "Bk": 97,
    "Cf": 98,
    "Es": 99,
    "Fm": 100,
    "Md": 101,
    "No": 102,
    "Lr": 103,
    "Rf": 104,
    "Db": 105,
    "Sg": 106,
    "Bh": 107,
    "Hs": 108,
    "Mt": 109,
    "Ds": 110,
    "Rg": 111,
    "Cn": 112,
    "Nh": 113,
    "Fl": 114,
    "Mc": 115,
    "Lv": 116,
    "Ts": 117,
    "Og": 118,
    "Uue": 119,
}
ATOMIC_NUMBERS = {value: key for key, value in ATOMIC_NUMBERS.items()}


def read_FEFF_geometry(path):
    path = Path(path) / "feff.inp"

    with open(path, "r") as f:
        feff_lines = f.readlines()
    where_atoms = [
        ii for ii, xx in enumerate(feff_lines) if xx.startswith("ATOMS")
    ]
    assert len(where_atoms) == 1
    feff_lines = feff_lines[where_atoms[0] + 3 : -1]
    feff_lines = [xx.split() for xx in feff_lines]

    atoms = [xx[4] for xx in feff_lines]

    distances = [float(xx[5]) for xx in feff_lines]

    # All together
    # 6 is the default precision used in FEFF by PMG!!!
    tmp = [
        [round(distance, ROUNDING_PRECISON), atom]
        for atom, distance in zip(atoms, distances)
    ]
    tmp.sort(key=lambda xx: xx[0])

    return {
        "atoms": [xx[1] for xx in tmp[1:]],
        "distances": np.array([xx[0] for xx in tmp[1:]]),
    }


def read_VASP_geometry(path, neighbor_radius=10.0):
    path = Path(path) / "POSCAR"

    # VASP POSCAR files are easy, only need data after line 8
    structure = IStructure.from_file(path)
    # vasp_coordinates = np.array([site.frac_coords for site in vasp_structure])
    # vasp_atoms = [site.specie.symbol for site in vasp_structure]

    neigh = structure.get_neighbors(structure[0], r=neighbor_radius)
    tmp = [
        [round(xx.nn_distance, ROUNDING_PRECISON), str(xx.specie)]
        for xx in neigh
    ]
    tmp.sort(key=lambda xx: xx[0])

    return {
        "atoms": [xx[1] for xx in tmp],
        "distances": np.array([xx[0] for xx in tmp]),
    }


def read_OCEAN_geometry(path, neighbor_radius=10.0):
    path = Path(path) / "ocean.in"

    absorber = str(path.parts[-2])

    with open(path, "r") as f:
        ocean_lines = f.readlines()
    ocean_lines = [xx.strip() for xx in ocean_lines]

    # Parse the ocean input file, it's kinda annoying...
    ocean_lines = (
        " ".join(ocean_lines).replace("{", "").replace("}", "|")[:-1].strip()
    )
    ocean_lines = ocean_lines.split("|")
    ocean_lines = [xx.strip() for xx in ocean_lines]

    geometry = None
    mat = None
    znucl = None
    typat = None

    for line in ocean_lines:
        split = line.split()
        key = split[0]

        if key == "xred":
            geometry = np.array(split[1:], dtype=float).reshape(-1, 3)
        elif key == "rprim":
            mat = np.array(split[1:], dtype=float).reshape(3, 3)
        elif key == "znucl":
            znucl = np.array(split[1:], dtype=int)
        elif key == "typat":
            typat = np.array(split[1:], dtype=int)

    assert geometry is not None
    assert mat is not None
    assert znucl is not None
    assert typat is not None

    atoms = [ATOMIC_NUMBERS[znucl[xx - 1]] for xx in typat]

    structure = IStructure(lattice=mat, species=atoms, coords=geometry)

    # Need to deal with this. OCEAN's indices of relevance are not necessarily 0
    absorbing_sites = [ii for ii, atom in enumerate(atoms) if atom == absorber]
    return_list = []
    for site in absorbing_sites:
        neigh = structure.get_neighbors(structure[site], r=neighbor_radius)
        tmp = [
            [round(xx.nn_distance, ROUNDING_PRECISON), str(xx.specie)]
            for xx in neigh
        ]
        tmp.sort(key=lambda xx: xx[0])
        return_list.append(
            {
                "atoms": [xx[1] for xx in tmp],
                "distances": np.array([xx[0] for xx in tmp]),
            }
        )
    return return_list


def read_XSpectra_geometry(path, neighbor_radius=10.0):
    path = Path(path) / "es.in"

    with open(path, "r") as f:
        xspectra_lines = f.readlines()
    xspectra_lines = [xx.strip() for xx in xspectra_lines]

    # Parse the Xspectra lines... also a pain, just like OCEAN
    atoms = []
    geom = []
    mat = None
    absorber = None
    ii = 0
    while ii < len(xspectra_lines):
        line = xspectra_lines[ii]
        if "ATOMIC_POSITIONS" in line:
            ii += 1
            jj = 0
            while True:
                line = xspectra_lines[ii]
                line_split = line.split()
                atom = line_split[0]
                position = line_split[1:]
                if "+" in atom:
                    absorber = jj
                    atom = atom[:-1]  # Remove the "+"
                if atom not in ATOMIC_NUMBERS.values():
                    break
                atoms.append(atom)
                geom.append(position)

                ii += 1
                jj += 1

        if "CELL_PARAMETERS" in line:
            mat = np.array(
                [xx.split() for xx in xspectra_lines[ii + 1 : ii + 4]],
                dtype=float,
            )
            break

        ii += 1

    assert absorber is not None

    geom = np.array(geom, dtype=float)

    structure = IStructure(lattice=mat, species=atoms, coords=geom)

    neigh = structure.get_neighbors(structure[absorber], r=neighbor_radius)
    tmp = [
        [round(xx.nn_distance, ROUNDING_PRECISON), str(xx.specie)]
        for xx in neigh
    ]
    tmp.sort(key=lambda xx: xx[0])

    return {
        "atoms": [xx[1] for xx in tmp],
        "distances": np.array([xx[0] for xx in tmp]),
    }


def read_EXCITING_geometry(path, neighbor_radius=10.0):
    path = Path(path) / "input.xml"
    structure = ExcitingInput.from_file(path).structure

    # Get the absorbing index
    with open(path, "r") as f:
        exciting_lines = f.readlines()
    exciting_lines = [xx.strip() for xx in exciting_lines]
    absorbing_line = [xx for xx in exciting_lines if "xasatom" in xx]
    assert len(absorbing_line) == 1
    absorbing_line = absorbing_line[0].split()
    absorbing_line = [xx for xx in absorbing_line if "xasatom" in xx]
    assert len(absorbing_line) == 1
    absorbing_line = absorbing_line[0]
    absorber = int(absorbing_line.split("=")[1].replace('"', "")) - 1

    neigh = structure.get_neighbors(structure[absorber], r=neighbor_radius)
    tmp = [
        [round(xx.nn_distance, ROUNDING_PRECISON), str(xx.specie)]
        for xx in neigh
    ]
    tmp.sort(key=lambda xx: xx[0])

    return {
        "atoms": [xx[1] for xx in tmp],
        "distances": np.array([xx[0] for xx in tmp]),
    }


@lru_cache(maxsize=16)
def _read_OCEAN_geometry(path, neighbor_radius):
    return read_OCEAN_geometry(path, neighbor_radius=neighbor_radius)


def _check_distances(
    d1,
    d2,
    d3,
    d4,
    a1,
    a2,
    a3,
    a4,
    path_FEFF,
    path,
    neighbor_radius,
    data_VASP,
    N,
):
    issues = []

    if not np.allclose(d1, d2, rtol=1e-5, atol=ATOL):
        issues.append(
            f"\nFEFF/VASP\n{len(a1)} | {a1}\n{len(a2)} | {a2}\n{d1}\n{d2}"
        )

    if not np.allclose(d2, d4, rtol=1e-5, atol=ATOL):
        issues.append(
            f"\nVASP/EXCITING\n{len(a2)} | {a2}\n{len(a4)} | {a4}\n{d2}\n{d4}"
        )

    if not np.allclose(d1, d4, rtol=1e-5, atol=ATOL):
        issues.append(
            f"\nFEFF/EXCITING\n{len(a1)} | {a1}\n{len(a4)} | {a4}\n{d1}\n{d4}"
        )

    if not np.allclose(d2, d3, rtol=1e-5, atol=ATOL):
        issues.append(
            f"\nVASP/Xspectra\n{len(a2)} | {a2}\n{len(a3)} | {a3}\n{d2}\n{d3}"
        )

    if not np.allclose(d1, d3, rtol=1e-5, atol=ATOL):
        issues.append(
            f"\nFEFF/Xspectra\n{len(a1)} | {a1}\n{len(a3)} | {a3}\n{d1}\n{d3}"
        )

    absorber = str(path_FEFF.name).split("_")[1]
    path_OCEAN = Path(path) / "OCEAN" / absorber

    all_data_OCEAN = _read_OCEAN_geometry(
        path_OCEAN, neighbor_radius=neighbor_radius
    )

    ocean_checks = []
    for data_OCEAN in all_data_OCEAN:
        if not (data_OCEAN["atoms"][:N] == data_VASP["atoms"][:N]):
            ocean_checks.append(False)
        else:
            ocean_checks.append(True)

        if not np.allclose(
            data_OCEAN["distances"][:N],
            data_VASP["distances"][:N],
            rtol=1e-5,
            atol=ATOL,
        ):
            ocean_checks.append(False)
        else:
            ocean_checks.append(True)

    if not any(ocean_checks):
        # raise AssertionError("OCEAN problem")
        issues.append("OCEAN problem")

    return issues


def consistency_check(path, neighbor_radius=10.0):
    """Summary

    Parameters
    ----------
    path : os.PathLike
        A path to a particular materials directory.
    """

    atom_dirs_FEFF = sorted(list((Path(path) / "FEFF").iterdir()))
    atom_dirs_VASP = sorted(list((Path(path) / "VASP").iterdir()))
    atom_dirs_XSpectra = sorted(list((Path(path) / "XSpectra").iterdir()))
    atom_dirs_EXCITING = sorted(list((Path(path) / "EXCITING").iterdir()))
    unit_cell_path = Path(path) / "POSCAR"

    # Check that the unit cell corresponds with the site indexes in the
    # VASP directory
    _sites = sorted([int(str(dd.name).split("_")[0]) for dd in atom_dirs_VASP])
    _structure = IStructure.from_file(unit_cell_path)
    inequivalent_sites = (
        SpacegroupAnalyzer(_structure)
        .get_symmetrized_structure()
        .equivalent_indices
    )
    for ineq in inequivalent_sites:
        L = len(ineq)
        L2 = len(list(set(ineq) - set(_sites)))
        assert L - L2 == 1
    # Done!

    l1 = [xx.name for xx in atom_dirs_FEFF]
    l2 = [xx.name for xx in atom_dirs_VASP]
    l3 = [xx.name for xx in atom_dirs_XSpectra]
    l4 = [xx.name for xx in atom_dirs_EXCITING]
    if not l1 == l2 == l4:
        raise AssertionError(f"\n{l1}\n{l2}\n{l4}")

    # Note that XSpectra comes with quite a few extra files...
    remaining_length = len(l3) - len(l2)
    if remaining_length != len(set(l3) - set(l2)):
        raise AssertionError(f"\n{l2}\n{l3}")

    for path_FEFF, path_VASP, path_XSpectra, path_EXCITING in zip(
        atom_dirs_FEFF, atom_dirs_VASP, atom_dirs_XSpectra, atom_dirs_EXCITING
    ):
        data_FEFF = read_FEFF_geometry(path_FEFF)
        data_VASP = read_VASP_geometry(
            path_VASP, neighbor_radius=neighbor_radius
        )
        data_XSpectra = read_XSpectra_geometry(
            path_XSpectra, neighbor_radius=neighbor_radius
        )
        data_EXCITING = read_EXCITING_geometry(
            path_EXCITING, neighbor_radius=neighbor_radius
        )

        a1 = data_FEFF["atoms"]
        N = len(a1)
        a2 = data_VASP["atoms"][:N]
        a3 = data_XSpectra["atoms"][:N]
        a4 = data_EXCITING["atoms"][:N]
        d1 = data_FEFF["distances"][:N]
        d2 = data_VASP["distances"][:N]
        d3 = data_XSpectra["distances"][:N]
        d4 = data_EXCITING["distances"][:N]

        # # Distances check
        issues = _check_distances(
            d1,
            d2,
            d3,
            d4,
            a1,
            a2,
            a3,
            a4,
            path_FEFF,
            path,
            neighbor_radius,
            data_VASP,
            N,
        )

        # TODO, if the atom types don't match, check the atom types in each
        # shell.
        # If the atom types don't line up, we definitely need to check
        # the distances
        # if not a1 == a2 == a3 == a4:
        #     issues.append(
        #         f"\nAtoms don't
        # match:\n{a1}\n{a2}\n{a3}\n{a4}\n{d1}\n{d2}\n{d3}\n{d4}"
        #     )

        if issues != []:
            raise AssertionError("\n".join(issues))

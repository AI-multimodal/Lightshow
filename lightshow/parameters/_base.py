from abc import ABC, abstractmethod
import numpy as np


def _get_k_mesh(structure, cutoff=32.0, max_radii=50.0):

    klist = dict()
    rlatt = np.array(structure.lattice.reciprocal_lattice.abc)

    # TODO: why 10, 0.2?
    for xx in np.arange(0, 10, 0.2):
        div = np.floor(xx * rlatt) + 1
        divlatt = 2.0 * np.pi / rlatt * div

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

    return k


class _BaseParameters(ABC):
    @abstractmethod
    def write(self, target_directory, kwargs):
        ...

    @property
    def name(self):
        return self._name

    @staticmethod
    def _getCondBands(volume, eRange):
        """Return a guess at the number of conduction bands that a given
        unit-cell volume needs to cover a given energy range (in Ryd).

        Parameters
        ----------
        volume : int
            The volume of the structure
        eRange : float
            The desired energy range

        Returns
        -------
        int :
            A number for the esimation number of bands
        """
        return round(0.256 * volume * (eRange ** (3 / 2)))

    @staticmethod
    def _getKmesh(structure, cutoff=32.0, max_radii=50.0):
        """Customizes the kmesh. For a kmesh sampling, e.g. [m, n, p],
        of a crystal "cell" is equivalent to generating a supercell with
        [m, n, p] the crystal cell. The corresponding radius is the largest
        radius of the sphere that can fit into this supercell. The radius can
        also be regarded as the inverse kmesh density. Along each reciprocal
        lattice, the kmesh densities are not guaranteed to be the same. The
        smallest kmesh density is chosen. This function uses the effective
        radius as a controlling factor to
        determining the kemsh.

        Parameters
        ----------
        supercell : pymatgen.core.structure.Structure
            Previously generated supercell (or standard unit cell).
        cutoff : float, optional
            Cutoff radius for constructing the kmesh. It will loop a look-up
            table for the effective radius (controlled by max_radii). The kmesh
            with radius right above the cutoff will be chosen. Default is 32
            Angstroms (60 Bohr).
        max_radii : float, optional
            Maximum radius used for constructing the lookup table.

        Returns
        -------
        tuple
        """

        return _get_k_mesh(structure, cutoff, max_radii)

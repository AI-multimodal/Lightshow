from abc import ABC, abstractmethod
from ase.units import Bohr

from monty.json import MSONable
import numpy as np


class _BaseKpointsMethod(ABC):
    @abstractmethod
    def __call__(self, structure):
        """Base abstraction for the various ways to get the kpoints for any
        simulation. This method should take the structure as input, and return
        a tuple containing the desired kmesh.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure

        Returns
        -------
        tuple
        """

        ...


class Constant(MSONable, _BaseKpointsMethod):
    def __init__(self, kmesh):
        self._kmesh = kmesh

    def __call__(self, structure):
        return tuple(self._kmesh)


class GenericEstimatorKpoints(MSONable, _BaseKpointsMethod):
    """For a kmesh sampling, e.g. [m, n, p], of a crystal "cell" is equivalent
    to generating a supercell with [m, n, p] the crystal cell. The
    corresponding radius is the largest radius of the sphere that can fit into
    this supercell. The radius can also be regarded as the inverse kmesh
    density. Along each reciprocal lattice, the kmesh densities are not
    guaranteed to be the same. The smallest kmesh density is chosen. This
    function uses the effective radius as a controlling factor to determining
    the kemsh.

    Parameters
    ----------
    cutoff : float, optional
        Cutoff radius for constructing the kmesh. It will loop a look-up
        table for the effective radius (controlled by max_radii). The kmesh
        with radius right above the cutoff will be chosen. Default is 32
        Angstroms (60 Bohr).
    max_radii : float, optional
        Maximum radius used for constructing the lookup table.
    """

    def __init__(self, cutoff=60.0, max_radii=80.0):
        self._cutoff = cutoff * Bohr
        self._max_radii = max_radii * Bohr

    def __call__(self, structure):
        klist = dict()
        rlatt = np.array(structure.lattice.reciprocal_lattice.abc)

        # TODO: why 10, 0.2?
        for xx in np.arange(0, 10, 0.2):
            div = np.floor(xx * rlatt) + 1
            divlatt = 2.0 * np.pi / rlatt * div

            radi = min(divlatt)

            if radi > self._max_radii:
                break
            else:
                div = tuple(div.astype(int))
                if div not in klist:
                    klist[div] = radi

        k = (1, 1, 1)
        for key, value in klist.items():
            if value > self._cutoff:
                k = key
                break

        return k

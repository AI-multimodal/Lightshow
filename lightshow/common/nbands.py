from abc import ABC, abstractmethod

from monty.json import MSONable


class _BaseNbandsMethod(ABC):
    @abstractmethod
    def __call__(self, structure):
        """Base abstraction for the various ways to get the number of bands
        used in a variety of alculations.

        Parameters
        ----------
        structure : pymatgen.core.structure.Structure

        Returns
        -------
        int
        """

        ...


class Constant(MSONable, _BaseNbandsMethod):
    def __init__(self, nbands):
        self._nbands = nbands

    def __call__(self, structure):
        return int(self._nbands)


class UnitCellVolumeEstimate(MSONable, _BaseNbandsMethod):
    """Return a guess at the number of conduction bands that a given
    unit-cell volume needs to cover a given energy range.

    Parameters
    ----------
    e_range : float
        The energy range in eV.
    """

    def __init__(self, e_range=40.0):
        self._e_range = e_range

    def __call__(self, structure):
        return int(
            round(
                0.256
                * structure.lattice.volume
                * ((self._e_range / 13.6056980659) ** (3.0 / 2.0))
            )
        )


# class LightshowHEG(MSONable, _BaseNbandsMethod):

#     def __init__(self, e_range=40.0):
#         self._e_range = e_range

#     def __call__(self, structure, total_valence_electrons):
#         """Gets an approximation for the total number of bands to use in the
#         calculation based on the free electron gas model, such that a range of
#         ``eRange`` is covered in the spectrum. The method assumes all
#         electrons
#         are free electrons, including those in the valence bands. The total
#         number of electrons will then be calculated as
#         :math:`N = N_\\mathrm{v} + N_\\mathrm{c}`, where :math:`N_\\mathrm{v}`
#         is the total number of valence electrons and
#         :math:`N_\\mathrm{c}` is the number of conduction band electrons above
#         the Fermi energy as specified by the provided energy range. The final
#         expression (in atomic units) is

#         .. math::

#             E_\\mathrm{range}/c_1 = c_2 / (V/N_vc * 3/ (4 pi))^2/3
#             - c_2 / (V/N_v * 3 / (4 pi))^2/3 N_v

#         where :math:`c_1 = 27.2114` and :math:`c_2 = 1.841`,
#         is summed over all valence electrons obtained from POTCAR.

#         Parameters
#         ----------
#         structure : pymatgen.core.structure.Structure
#         eRange : float
#             The approximate energy range desired in the spectrum, starting at
#             the Fermi energy.

#         Returns
#         -------
#         int
#             Approximation for the number of bands to use in computing the
#             XANES
#             spectrum as computed using the (3D) homogeneous electron gas
#             approximation.
#         """

#         ANG2BOHR = 1.88873
#         volume_bohr = structure.volume * ANG2BOHR**3
#         Nv = self.get_total_valence_electrons(structure)

#         tmp1 = self._e_range / 27.2114 + 1.841 / (
#             volume_bohr / total_valence_electrons * 3.0 / 4.0 / np.pi
#         ) ** (2 / 3)
#         tmp2 = (1.841 / tmp1) ** (2 / 3) / 3.0 * 4.0 * np.pi
#         Nvc = volume_bohr / tmp2 * 0.65
#         return int(np.ceil(Nvc))

"""XAS broadening using manual or empirically optimized parameters"""


from scipy.stats import pearsonr as ps
from scipy import optimize
from scipy.stats import norm
from scipy.stats import cauchy
from scipy.special import voigt_profile as voigt

import numpy as np


def gauss_broaden(x, xin, yin, sigma, shift=0.0):
    """Gaussian broadening function

    Parameters
    ----------
    x : numpy.ndarray
        The output energy grid.
    xin : numpy.ndarray
        The input energy grid.
    yin : numpy.ndarray
        The input spectrum.
    sigma : float
        Gaussian broadening sigma.
    shift : float
        The energy shift.

    Returns
    -------
    numpy.ndarray
        Broadened spectrum
    """

    x1, x2 = np.meshgrid(x, xin + shift)

    return (
        np.dot(norm.pdf(x1, x2, sigma).T, yin) / len(xin) * (xin[-1] - xin[0])
    )


def lorentz_broaden(x, xin, yin, gamma, shift=0.0):
    """Lorentzian broadening function

    Parameters
    ----------
    x : numpy.ndarray
        The output energy grid.
    xin : numpy.ndarray
        The input energy grid.
    yin : numpy.ndarray
        The input spectrum.
    gamma : float
        Lorentzian broadening gamma.
    shift : float
        The energy shift.

    Returns
    -------
    numpy.ndarray
        Broadened spectrum.
    """

    x1, x2 = np.meshgrid(x, xin + shift)

    return (
        np.dot(cauchy.pdf(x1, x2, gamma).T, yin) / len(xin) * (xin[-1] - xin[0])
    )


def voigt_broaden(x, xin, yin, sigma, gamma, shift=0.0):
    """Voigt broadening function

    Parameters
    ----------
    x : numpy.ndarray
        The output energy grid.
    xin : numpy.ndarray
        The input energy grid.
    yin : numpy.ndarray
        The input spectrum.
    sigma : float
        Voigt broadening sigma.
    gamma : float
        Voigt broadening gamma.
    shift : float
        The energy shift.

    Returns
    -------
    numpy.ndarray
        Broadened spectrum.
    """

    x1, x2 = np.meshgrid(x, xin + shift)

    return (
        np.dot(voigt(x1 - x2, sigma, gamma / 2).T, yin)
        / len(xin)
        * (xin[-1] - xin[0])
    )


class SiteSpectum:
    """This class stores information about an XAS spectrum at each site,
    including the spectrum itself, the Fermi energy and site weight.

    Parameters
    ----------
    xin : numpy.ndarray
        The input energy grid.
    yin : numpy.ndarray
        The input spectrum.
    fermi : float
        The fermi energy of the spectrum.
    weight : int
        The weight of the site spectrum.
    """

    def __init__(
        self,
        xin,
        yin,
        fermi=0.0,
        weight=1,
    ):
        self._xin = xin
        self._yin = yin
        self._Fermi = fermi
        self._Weight = weight

    @property
    def spectrum(self):
        return self._xin, self._yin

    @property
    def weight(self):
        return self._Weight

    @property
    def fermi(self):
        return self._Fermi


class Broaden:
    """The class that store broadening parameters and provide functions to
    optimize broadening parameters empirically.

    Parameters
    ----------
    sigma : float
        The gaussian broadening sigma.
    lorentz_divider : float
        The parameter to determine the energy dependent gamma with equation:
        :math:`\\Gamma\\left(E\right)=\\Gamma_0+(E-E_f)/l`,
        where :math:`l` is the lorentz_divider, and the :math:`\\Gamma_0` is
        determined by Corehole life time.
    CHlifetime : float
        Core hole life time, e.g. Zn K-edge: 1.67, Ti K-edge: 0.89.
    shift : float
        The energy shift.
    """

    def __init__(
        self, sigma=0.0, lorentz_divider=1e5, CHlifetime=1.0, shift=0.0
    ):
        self._sigma = sigma
        self._lorentz_divider = lorentz_divider
        self._CHlifetime = CHlifetime
        self._shift = shift

    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, sigma):
        self._sigma = sigma

    @property
    def lorentz_divider(self):
        return self._lorentz_divider

    @lorentz_divider.setter
    def lorentz_divider(self, lorentz_divider):
        self._lorentz_divider = lorentz_divider

    @property
    def shift(self):
        return self._shift

    @shift.setter
    def shift(self, shift):
        self._shift = shift

    # Do energy dependent voigt function broadening
    @staticmethod
    def _ene_dep_broaden(x, xin, yin, sigma, ld, lt, ef):
        x1, x2 = np.meshgrid(x, xin)

        gam = lt + np.heaviside(x2 - ef, 1) * (x2 - ef) / ld

        return (
            np.dot(voigt(x1 - x2, sigma, gam / 2).T, yin)
            / len(xin)
            * (xin[-1] - xin[0])
        )

    # Broaden spectrum from each site and combine them
    def _spectrum_broaden(
        self, x, site_spectra, sigma, ld, shift, cross_sec=True
    ):
        yout = np.zeros_like(x)

        weight = 0

        for site_spectrum in site_spectra:
            xin, yin = site_spectrum.spectrum

            xin = xin + shift

            if cross_sec:
                yin = yin * xin

            yout = (
                yout * weight
                + self._ene_dep_broaden(
                    x,
                    xin,
                    yin,
                    sigma,
                    ld,
                    self._CHlifetime,
                    site_spectrum.fermi,
                )
                * site_spectrum.weight
            )

            weight = weight + site_spectrum.weight

        return yout / weight

    # Give the pearson correlation score between broadened simulation and
    # experimental spectrum
    def _broaden_score(
        self,
        broaden_paras,
        x,
        site_spectra,
        exp,
        dmu=True,
        cross_sec=True,
        *args,
    ):
        spec = self._spectrum_broaden(
            x,
            site_spectra,
            broaden_paras[0],
            broaden_paras[1],
            broaden_paras[2],
            cross_sec,
        )

        if dmu:
            exp = np.diff(exp) / np.diff(x)

            spec = np.diff(spec) / np.diff(x)

        return -ps(spec, exp)[0]

    def broaden(self, x, site_spectra, volume=1.0, cross_sec=True):
        """Process broadening with the energy dependent voigt function and the
        stored parameters.

        Parameters
        ----------
        x : numpy.ndarray
            The output energy grid.
        site_spectra : list
            The list of site_spectrum objects
        volume : float
            The factor that tune the intensity of broadened spectrum.
        cross_sec : bool
            If True, the spectrum will be converted to be absorption cross
            section.

        Returns
        -------
        numpy.ndarray
        """

        return (
            self._spectrum_broaden(
                x,
                site_spectra,
                self._sigma,
                self._lorentz_divider,
                self._shift,
                cross_sec,
            )
            * volume
        )

    def paras_optimize(
        self,
        x,
        site_spectra,
        exp,
        bounds=None,
        opt_shift=False,
        dmu=True,
        cross_sec=True,
    ):
        """Optimize the broadening parameter with simplicial homology global
        optimization in scipy package.

        Parameters
        ----------
        x : numpy.ndarray
            The output energy grid.
        site_spectra : list
            The list of site_spectrum objects.
        exp : numpy.ndarray
            The experimental spectrum.
        bounds : list
            The bounds of the broadening parameters.
        opt_shift : bool
            If True, the method will optimize energy shift value.
        dmu : bool
            If True, the method will compare the first derivatives of
            experimental spectrum and simulation.
        cross_sec : bool
            If true, the spectrum will be converted to be absorption cross
            section.

        Returns
        -------
        numpy.ndarray
        """

        if opt_shift is False:
            xi = 10e-3

            bounds[2] = (self._shift - xi, self._shift + xi)

        result = optimize.shgo(
            self._broaden_score,
            bounds,
            args=(x, site_spectra, exp, dmu, cross_sec),
        )

        self._sigma = result.x[0]
        self._lorentz_divider = result.x[1]

        if opt_shift is True:
            self._shift = result.x[2]

        return result

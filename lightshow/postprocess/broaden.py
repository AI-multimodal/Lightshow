"""XAS broadening using manual or empirically optimized parameters"""

from monty.json import MSONable
import numpy as np
from scipy.stats import pearsonr
from scipy import optimize
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.stats import norm
from scipy.stats import cauchy
from scipy.special import voigt_profile as voigt


def gauss_broaden(xout, xin, yin, sigma):
    """Applies Gaussian broadening to the input data, with standard deviation
    sigma.

    .. math::

        G(x; \\sigma) = \\frac{1}{\\sqrt{2 \\pi \\sigma^2}}
        e^{-x^2 / 2\\sigma^2}

    Parameters
    ----------
    xout : numpy.ndarray
        The output x-axis grid.
    xin : numpy.ndarray
        The input x-axis grid.
    yin : numpy.ndarray
        The input data for f(x) to broaden. Generally a spectrum.
    sigma : float
        The standard deviation of the Gaussian function used for broadening.

    Returns
    -------
    numpy.ndarray
        Broadened spectrum
    """

    x1, x2 = np.meshgrid(xout, xin)
    dx = xin[-1] - xin[0]
    return np.dot(norm.pdf(x1, x2, sigma).T, yin) / len(xin) * dx


def lorentz_broaden(x, xin, yin, gamma):
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
        Lorentzian broadening Gamma or the full-width at half-maximum (FWHM)

    Returns
    -------
    numpy.ndarray
        Broadened spectrum.
    """

    x1, x2 = np.meshgrid(x, xin)
    dx = xin[-1] - xin[0]
    return np.dot(cauchy.pdf(x1, x2, gamma / 2.0).T, yin) / len(xin) * dx


def voigt_broaden(x, xin, yin, sigma, gamma):
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

    Returns
    -------
    numpy.ndarray
        Broadened spectrum.
    """

    x1, x2 = np.meshgrid(x, xin)
    dx = xin[-1] - xin[0]
    return np.dot(voigt(x1 - x2, sigma, gamma / 2.0).T, yin) / len(xin) * dx


def energy_dependent_voigt_broaden(x, xin, yin, sigma, alpha, lifetime, efermi):
    """See ``voigt_broaden``. Performs a similar broadening only with an
    energy-dependent Lorentzian gamma term.

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
    alpha : float
        Parameter alpha in the energy dependent broadening function
    lifetime : float
        Core-hole lifetime
    efermi : float
        Fermi energy

    Returns
    -------
    numpy.ndarray
    """

    x1, x2 = np.meshgrid(x, xin)
    gamma = lifetime + np.heaviside(x2 - efermi, 1) * (x2 - efermi) * alpha
    dx = xin[-1] - xin[0]
    return np.dot(voigt(x1 - x2, sigma, gamma / 2.0).T, yin) / len(xin) * dx


def get_shift_by_maxima(exp, theory):
    exp_x = exp[:, 0].copy()
    exp_y = exp[:, 1].copy()
    theory_x = theory[:, 0].copy()
    theory_y = theory[:, 1].copy()

    # Start with an initial guess: we shift the theory_x location by the
    # difference between the location of the max experimental peak and max
    # theory peak
    exp_max_x = exp_x[np.argmax(exp_y)]
    thr_max_x = theory_x[np.argmax(theory_y)]
    return thr_max_x - exp_max_x


def align_spectra_by_maxima(exp, theory):
    """Performs a simple alignment of the theoretical spectrum to the
    experiment by using the peak maxima. Returns the new x-axis.

    Parameters
    ----------
    exp : numpy.ndarray
    theory : numpy.ndarray

    Returns
    -------
    numpy.ndarray
        The aligned x-axis grid.
    """

    return theory[:, 0].copy() - get_shift_by_maxima(exp, theory)


class SiteSpectrum(MSONable):
    """This class stores information about an XAS spectrum at each site,
    including the spectrum itself, the Fermi energy and site weight.

    Parameters
    ----------
    spectrum : numpy.ndarray
        The N x 2 array corresponding to the spectrum.
    e_fermi : float
        The e_fermi energy of the spectrum.
    weight : int
        The weight of the site spectrum.
    """

    @property
    def original_spectrum(self):
        return self._original_spectrum

    @property
    def spectrum(self):
        return self._spectrum

    @property
    def weight(self):
        return self._weight

    @property
    def e_fermi(self):
        return self._e_fermi

    def __init__(
        self,
        spectrum,
        original_spectrum=None,
        e_fermi=None,
        weight=1.0,
        broadening_kwargs=None,
    ):
        self._spectrum = spectrum
        self._original_spectrum = original_spectrum
        if self._original_spectrum is None:
            self._original_spectrum = self._spectrum.copy()
        self._e_fermi = e_fermi
        self._weight = weight
        self._broadening_kwargs = broadening_kwargs

    def scale_max_to_1_(self):
        self._spectrum[:, 1] /= self._spectrum[:, 1].max()

    def align_to_experiment_(self, exp):
        """Modifies the spectrum in-place to align it to the provided
        experimental result. Note that this destroys the original information
        about the spectrum, and assumes that the original energy grid of the
        theory data is unusable.

        Parameters
        ----------
        exp : numpy.ndarray
            The experimental data to align to.
        """

        self._spectrum[:, 0] = align_spectra_by_maxima(exp, self._spectrum)

    def broaden(
        self, xout, method="Gaussian", *broadening_args, **broadening_kwargs
    ):
        """Returns a broadened version of the spectrum.

        Parameters
        ----------
        xout : TYPE
            Description
        method : str, optional
            The method to use. Should either match the function signature or
            be "Gaussian" or "Lorentzian".
        *broadening_args
            Description
        **broadening_kwargs
            Extra keyword arguments to provide for the broadening. Some of
            these will be provided. Note for "energy_dependent_voigt_broaden"
            e_fermi is provided as an attribute.

        Returns
        -------
        numpy.ndarray
        """

        xin = self._spectrum[:, 0].copy()  # + shift
        yin = self._spectrum[:, 1].copy()

        if method == "Gaussian" or method == "gauss_broaden":
            res = gauss_broaden(
                xout, xin, yin, *broadening_args, **broadening_kwargs
            )
        elif method == "Lorentzian" or method == "lorentz_broaden":
            res = lorentz_broaden(
                xout, xin, yin, *broadening_args, **broadening_kwargs
            )
        elif method == "Voigt" or method == "voigt_broaden":
            res = voigt_broaden(
                xout, xin, yin, *broadening_args, **broadening_kwargs
            )
        elif method == "energy_dependent_voigt_broaden":
            res = energy_dependent_voigt_broaden(
                xout,
                xin,
                yin,
                *broadening_args,
                efermi=self._e_fermi,
                **broadening_kwargs,
            )
        else:
            raise ValueError(f"Uknown method {method}")
        self._broadening_kwargs = broadening_kwargs
        return res

    def optimize_parameters(
        self,
        exp,
        bounds,
        method="energy_dependent_voigt_broaden",
    ):
        """Optimizes the broadening parameters for a single spectrum. A
        reference spectrum (``exp``) is provided as the ground truth.

        Parameters
        ----------
        exp : numpy.ndarray
            Description
        bounds : None, optional
            Note that shift is always the first parameter.
        method : str, optional
            Description

        """

        xin_theory = self._spectrum[:, 0].copy()

        # Bounds[0][0] is negative
        assert bounds[0][0] <= 0.0
        lower_shift = exp[0, 0] + bounds[0][0]
        upper_shift = exp[-1, 0] + bounds[0][1]

        self.align_to_experiment_(exp)  # makes the shift optimization easier

        where_lower = np.argmin((xin_theory - lower_shift) ** 2)
        where_upper = np.argmin((xin_theory - upper_shift) ** 2)

        xin_theory = xin_theory[where_lower:where_upper]

        def objective(params):
            shift = params[0]
            broadened_spectrum = self.broaden(xin_theory, method, *params[1:])
            ius = InterpolatedUnivariateSpline(
                xin_theory + shift, broadened_spectrum, k=3
            )
            return -pearsonr(exp[:, 1], ius(exp[:, 0]))[0]

        result = optimize.shgo(objective, bounds)
        return result

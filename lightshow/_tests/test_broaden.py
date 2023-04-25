import pytest

import numpy as np

from lightshow.postprocess.broaden import SiteSpectrum


@pytest.fixture
def anatase_exp(spectra_file_directory):
    # Anatase TiO2 experimental data is derived from Yan, D.,et al.
    # Nano letters, 19(6), 3457-3463.
    return np.loadtxt(spectra_file_directory / "anatase_exp.txt")


@pytest.fixture
def anatase_VASP_theory(spectra_file_directory):
    # Anatase TiO2 theoretical simulation is calculated by VASP
    unbroaden_spect = np.loadtxt(
        spectra_file_directory / "anatase_theory_VASP.txt"
    )

    return unbroaden_spect[30000:35000, :]


def test_broaden(anatase_VASP_theory, anatase_exp, sample_file_directory):
    """Makes sure the broadening code at least works!"""

    # Some information we'll need for the energy-dependent Voigt broadening
    TiK_core_state = -4864.0371
    TiK_core_hole_lifetime = 0.89
    ana_core_hole_fermi_energy = 5.1302
    e_fermi = ana_core_hole_fermi_energy - TiK_core_state

    # Output energy grid
    x = np.linspace(4950, 5010, 1000, endpoint=True)

    vasp_spectrum = SiteSpectrum(anatase_VASP_theory.copy(), e_fermi=e_fermi)
    vasp_spectrum.align_to_experiment_(anatase_exp)
    vasp_spectrum.scale_max_to_1_()

    vasp_spectrum.broaden(x, method="Gaussian", sigma=1.0)

    # ... using Lorentzian broadening
    vasp_spectrum.broaden(x, method="Lorentzian", gamma=1.0)

    # ... using Voigt broadening
    vasp_spectrum.broaden(x, method="Voigt", sigma=0.5, gamma=0.5)

    # ... using energy-dependent Voigt broadening
    vasp_spectrum.broaden(
        x,
        method="energy_dependent_voigt_broaden",
        sigma=0.2,
        alpha=0.025,
        lifetime=TiK_core_hole_lifetime,
    )

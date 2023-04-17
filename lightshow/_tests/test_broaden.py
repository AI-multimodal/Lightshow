import pytest

import numpy as np

from lightshow.postprocess import broaden


@pytest.fixture
def anatase_exp(spectra_file_directory):
    # Anatase TiO2 experimental data is derived from Yan, D.,et al.
    # Nano letters, 19(6), 3457-3463.
    exp = np.loadtxt(spectra_file_directory / "anatase_exp.txt")
    return exp.T


@pytest.fixture
def anatase_theory(spectra_file_directory):
    # Anatase TiO2 theoretical simulation is calculated by VASP
    unbroaden_spect = np.loadtxt(spectra_file_directory / "anatase_theory.txt")

    # choosing smaller range can make the optimization processs faster
    xin = unbroaden_spect.T[0][20000:35000]
    yin = unbroaden_spect.T[1][20000:35000]
    return [xin, yin]


def test_broaden(anatase_theory, sample_file_directory):
    # Output energy grid
    x = np.linspace(4950, 5010, 1000, endpoint=True)

    # Create site spectrum object for broadening
    ana = broaden.SiteSpectum(anatase_theory[0], anatase_theory[1], weight=1.0)

    spect = broaden.gauss_broaden(
        x, ana.spectrum[0], ana.spectrum[1], sigma=1.0, shift=100
    )
    assert (
        max(
            abs(
                spect
                - np.loadtxt(sample_file_directory / "gauss_broad.txt").T[1]
            )
        )
        < 1e-10
    )

    spect = broaden.lorentz_broaden(
        x, ana.spectrum[0], ana.spectrum[1], gamma=1.0, shift=100
    )
    assert (
        max(
            abs(
                spect
                - np.loadtxt(sample_file_directory / "lorentz_broad.txt").T[1]
            )
        )
        < 1e-10
    )

    spect = broaden.voigt_broaden(
        x, ana.spectrum[0], ana.spectrum[1], sigma=0.5, gamma=0.5, shift=100
    )
    assert (
        max(
            abs(
                spect
                - np.loadtxt(sample_file_directory / "voigt_broad.txt").T[1]
            )
        )
        < 1e-10
    )


def test_paras_optimize(anatase_theory, anatase_exp, sample_file_directory):
    # Collect XAS sumulation information
    TiK_core_state = -4864.0371
    TiK_core_hole_lifetime = 0.89

    Ana_core_hole_fermi_energy = 5.1302

    Fermi = Ana_core_hole_fermi_energy - TiK_core_state

    # Output energy grid
    x = np.linspace(4950, 5010, 1000, endpoint=True)

    # Enter the information to site spectrum object
    ana = broaden.SiteSpectum(
        anatase_theory[0], anatase_theory[1], fermi=Fermi, weight=1.0
    )

    brd = broaden.Broaden(
        sigma=0.2,
        lorentz_divider=40,
        CHlifetime=TiK_core_hole_lifetime,
        shift=100,
    )

    spect = brd.broaden(x, [ana])
    assert (
        max(
            abs(
                spect
                - np.loadtxt(sample_file_directory / "ene_dep_broad.txt").T[1]
            )
        )
        < 1e-10
    )

    bds = [(0.1, 2), (1, 50), (95, 105)]

    mes = brd.paras_optimize(
        anatase_exp[0],
        [ana],
        anatase_exp[1],
        bounds=bds,
        dmu=True,
        opt_shift=False,
    )
    assert mes.message == "Optimization terminated successfully."

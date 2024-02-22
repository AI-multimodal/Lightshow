# This is a simple tool to compare to plots
# Chuntian Cao developed based on Fanchen Meng's 2022 code


import numpy as np
from scipy import interpolate
from scipy.stats import pearsonr, spearmanr


def compare_between_spectra(
    spectrum1, spectrum2, erange=35, accuracy=0.01, method="coss"
):
    """Atomatic align the spectra and calculate the correlation coefficients.

    Parameters
    ----------
    spectrum1, spectrum2 : 2d-array
        Two-column arrays of energy vs. intensity XAS data.
    method : {'coss', 'pearson', 'spearman'}
        The correlation metric for spectra comparison.
        Empirically 'coss' works well.
    erange : float, default=35
        Energy range for comparison. Unit: eV.
    accuracy : float, default=0.01
        Accuracy for spectra alignment. Unit: eV.

    Returns
    -------
    pearson, spearman, coss : float
        Correlation of two spectra after alignment.
    shift : float
        Relative shift between the two spectra, sign is meaningful.
        Spectrum2 should be shifted to spectrum2+shift for alignment.

    """

    start1, end1 = truncate_spectrum(spectrum1, erange)
    plot1 = np.column_stack(
        (
            spectrum1[start1:end1, 0] - spectrum1[start1][0],
            spectrum1[start1:end1, 1],
        )
    )
    start2, end2 = truncate_spectrum(spectrum2, erange)
    plot2 = np.column_stack(
        (
            spectrum2[start2:end2, 0] - spectrum2[start2][0],
            spectrum2[start2:end2, 1],
        )
    )

    shift_prior = spectrum1[start1, 0] - spectrum2[start2, 0]
    _, shift = max_corr(plot1, plot2, step=accuracy, method=method)
    pears, spear, coss = spectra_corr(
        plot1, plot2, omega=shift, verbose=True, method="all"
    )
    shift += shift_prior

    return pears, spear, coss, shift


def truncate_spectrum(spectrum, erange=35, threshold=0.02):
    """Truncate XAS spectrum to desired energy range.

    Parameters
    ----------
    spectrum : 2d-array
        Column stacked spectrum, energy vs. intensity.
    erange : float, default=35
        Truncation energy range in eV.
    threshold : float, default=0.02
        Start truncation at threshold * maximum intensity.

    Returns
    -------
    start, end : int
        Indices of truncated spectrum in the input spectrum.

    """
    x = spectrum[:, 0]
    y = spectrum[:, 1] / np.max(spectrum[:, 1])

    logic = y > threshold
    seq = x == x[logic][0]
    start = seq.argmax()

    logic = x < x[start] + erange
    seq = x == x[logic][-1]
    end = seq.argmax()

    return start, end


def cos_similar(v1, v2):
    """Calculates the cosine similarity between two vectors.

    Parameters
    ----------
    v1, v2 : 1d-array

    Returns
    -------
    cosSimilarity : float

    """
    norm1 = np.sqrt(np.dot(v1, v1))
    norm2 = np.sqrt(np.dot(v2, v2))
    cosSimilarity = np.dot(v1, v2) / (norm1 * norm2)
    return cosSimilarity


def spectra_corr(
    spectrum1, spectrum2, omega=0, grid=None, verbose=True, method="all"
):
    """Calculates pearson, spearman, and cosine correlation of two spectra.

    Parameters
    ----------
    spectrum1, spectrum2 : 2d-array
        Two-column arrays of energy vs. intensity XAS data.
    omega : float
        Shift between two spectra. spectrum2 shifted to spectrum2 + omega.
    grid : 1d-array
        Common grid for interpolation.
    method : {'all', 'pearson', 'spearman', 'coss'}

    Returns
    -------
    correlation : list of float
        Pearson, spearman, or cosine similarity.
        If method == 'all', all three correlations are calculated.

    """
    if grid is None:
        grid = np.linspace(
            max(spectrum1[0, 0], spectrum2[0, 0] + omega),
            min(spectrum1[-1, 0], spectrum2[-1, 0] + omega),
            300,
        )
    interp1 = interpolate.interp1d(
        spectrum1[:, 0],
        spectrum1[:, 1],
        assume_sorted=False,
        kind="cubic",
        bounds_error=False,
    )
    interp2 = interpolate.interp1d(
        spectrum2[:, 0] + omega,
        spectrum2[:, 1],
        assume_sorted=False,
        kind="cubic",
        bounds_error=False,
    )
    curve1 = interp1(grid)
    curve2 = interp2(grid)
    indices = ~(np.isnan(curve1) | np.isnan(curve2))

    pearson = np.NaN
    spearman = np.NaN
    coss = np.NaN
    if "pear" in method:
        pearson = pearsonr(curve1[indices], curve2[indices])[0]
    if "spear" in method:
        spearman = spearmanr(curve1[indices], curve2[indices])[0]
    if "cos" in method:
        coss = cos_similar(curve1[indices], curve2[indices])
    if method == "all":
        pearson = pearsonr(curve1[indices], curve2[indices])[0]
        spearman = spearmanr(curve1[indices], curve2[indices])[0]
        coss = cos_similar(curve1[indices], curve2[indices])
    width = 0.5 * min(
        spectrum1[-1, 0] - spectrum1[0, 0], spectrum2[-1, 0] - spectrum2[0, 0]
    )
    # require 50% overlap

    if grid[indices][-1] - grid[indices][0] < width:
        decay = 0.9 ** (width / (grid[indices][-1] - grid[indices][0]))
        if verbose:
            print(
                "Overlap less than 50%%. Similarity values decayed by %0.4f"
                % decay
            )
        pearson *= decay
        spearman *= decay
        coss *= decay
    correlation = [ii for ii in [pearson, spearman, coss] if not np.isnan(ii)]
    return correlation


def max_corr(
    spectrum1,
    spectrum2,
    start=12,
    stop=-12,
    step=0.01,
    grid=None,
    method="coss",
):
    """Calculate the correlation between two spectra,
        and the amout of shift to obtain maximum correlation.

    Parameters
    ----------
    spectrum1, spectrum2 : 2d-array
        Two-column arrays of energy vs. intensity XAS.
    start, stop, step : float
        Shift of spectrum2 ranges from start to stop with stepsize=step.
    grid : 1d-array
        Common grid for interpolation.
    method : {'coss', 'pearson', 'spearman'}
        Empirically 'coss' (cosine similarity) works well.

    Returns
    -------
    correlation : dict
        Correlation values at each shift step.
    m_shift : float
        Shift value at which the correlation is max.

    """

    if start <= stop:
        print("WARNING: Start {} is larger than stop {}]".format(start, stop))
        exit()
    correlation = {}

    i = start
    while i > stop:
        correlation[i] = spectra_corr(
            spectrum1,
            spectrum2,
            omega=i,
            grid=grid,
            verbose=False,
            method=method,
        )[0]
        i -= step
    # find index at maximum correlation

    m = 0
    for i, j in correlation.items():
        if j > m:
            m = j
            m_shift = i
    # check if the gradient makes sense

    gplot1 = np.vstack(
        (
            spectrum1[:, 0],
            np.gradient(spectrum1[:, 1], spectrum1[1, 0] - spectrum1[0, 0]),
        )
    ).T
    gplot2 = np.vstack(
        (
            spectrum2[:, 0],
            np.gradient(spectrum2[:, 1], spectrum2[1, 0] - spectrum2[0, 0]),
        )
    ).T
    x1 = peak_loc(gplot1)
    x2 = peak_loc(gplot2)
    if abs(x1 - m_shift - x2) < 2:
        pass
    else:
        print(
            "XAS edge positions might not align. "
            "Better to plot and check the spectrum."
        )
    return correlation, m_shift


def peak_loc(plot):
    """Locate the peak positon of a spectrum.

    Parameters
    ----------
    plot : 2d-array

    Returns
    -------
    position of the peak

    """
    return plot[plot[:, 1].argmax(), 0]

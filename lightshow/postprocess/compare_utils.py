# This is a simple tool to compare to plots
# Chuntian Cao developed based on Fanchen Meng's 2022 code


import numpy as np
from scipy import interpolate
from scipy.stats import pearsonr, spearmanr


def compare_between_spectra(
    spectrum1, spectrum2, erange=35, accuracy=0.01, method="coss"
):
    """Atomatic align the spectra and calculate the spearman coefficient
    Comparison was done for epsilon instead of cross-section

    Parameters
    ----------
    spectrum1 and spectrum2 : two-column arrays of energy vs. intensity XAS
    method : 'pearson', 'spearman', or 'coss' (cosine similarity). 
        Empirically 'coss' works well

    Return
    ------
    pearson, spearman, cossine correlation : real
    shift : real; relative shift between the two spectra, sign is meaningful
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
    _, shift = maxCorr(plot1, plot2, step=accuracy, method=method)
    pears, spear, coss = spectraCorr(
        plot1, plot2, omega=shift, verbose=True, method="all"
    )
    shift += shift_prior

    return pears, spear, coss, shift


def truncate_spectrum(spectrum, erange=35, threshold=0.02):
    """
    Parameters
    ----------
    spectrum: column stacked spectrum, energy vs. intensity
    erange: truncation range
    threshold: start truncation at threshold * maximum intensity

    Return
    ------
    start, end : int; indices of truncated spectrum in the input spectrum
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


def CosSimilar(v1, v2):
    """Calculates the cosine similarity between two vectors

    Parameters
    ----------
    v1 : vector #1
    v2 : vector #2

    Returns
    -------
    cosine similarity : real; v1 * v2 / (norm(v1) * norm(v2))
    """
    norm1 = np.sqrt(np.dot(v1, v1))
    norm2 = np.sqrt(np.dot(v2, v2))
    cosSimilarity = np.dot(v1, v2) / (norm1 * norm2)
    return cosSimilarity


def spectraCorr(
    spectrum1, spectrum2, omega=0, GRID=None, verbose=True, method="all"
):
    """Calculates pearson, spearman, cossine correlation between two spectra

    Parameters
    ----------
    spectrum1 and spectrum2 : two-column arrays of energy vs. intensity XAS
    omega : shift between two spectra. spectrum2 shifted to spectrum2 + omega
    GRID : common grid for interpolation
    method : 'pearson', 'spearman', 'coss' (cosine similarity), or 'all'.

    Returns
    -------
    correlation: list; pearson, spearman, or cosine similarity. 
        If method == 'all', all three correlations are returned.
    """
    if GRID is None:
        GRID = np.linspace(
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
    curve1 = interp1(GRID)
    curve2 = interp2(GRID)
    indices = ~(np.isnan(curve1) | np.isnan(curve2))

    pearson = np.NaN
    spearman = np.NaN
    coss = np.NaN
    if "pear" in method:
        pearson = pearsonr(curve1[indices], curve2[indices])[0]
    if "spear" in method:
        spearman = spearmanr(curve1[indices], curve2[indices])[0]
    if "cos" in method:
        coss = CosSimilar(curve1[indices], curve2[indices])
    if method == "all":
        pearson = pearsonr(curve1[indices], curve2[indices])[0]
        spearman = spearmanr(curve1[indices], curve2[indices])[0]
        coss = CosSimilar(curve1[indices], curve2[indices])
    width = 0.5 * min(
        spectrum1[-1, 0] - spectrum1[0, 0], spectrum2[-1, 0] - spectrum2[0, 0]
    )
    # require 50% overlap

    if GRID[indices][-1] - GRID[indices][0] < width:
        decay = 0.9 ** (width / (GRID[indices][-1] - GRID[indices][0]))
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


def maxCorr(
    spectrum1,
    spectrum2,
    start=12,
    stop=-12,
    step=0.01,
    GRID=None,
    method="coss",
):
    """Calculate the correlation between two spectra, 
        and the amout of shift to obtain maximum correlation

    Parameters
    ----------
    spectrum1 and spectrum2 : two-column arrays of energy vs. intensity XAS
    start, stop, and step : shift of spectrum2 ranges from start to stop 
        with stepsize=step
    GRID : common grid for interpolation
    method : 'pearson', 'spearman', or 'coss' (cosine similarity). 
        Empirically 'coss' works well

    Returns
    -------
    correlation : dict; values at each shift step
    m_shift : the shift value at which the correlation is max

    """

    if start <= stop:
        print("WARNING: Start {} is larger than stop {}]".format(start, stop))
        exit()
    correlation = {}

    i = start
    while i > stop:
        correlation[i] = spectraCorr(
            spectrum1,
            spectrum2,
            omega=i,
            GRID=GRID,
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
            "XAS edge positions might not align. " \
            "Better to plot and check the spectrum. "
        )
    return correlation, m_shift


def peak_loc(plot):
    """locate the peak positon of a spectrum

    Parameters
    ----------
    plot : 2d-array

    Returns
    -------
    position of the peak
    """
    return plot[plot[:, 1].argmax(), 0]

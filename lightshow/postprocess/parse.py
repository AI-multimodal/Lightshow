import numpy as np
from scipy.constants import physical_constants
import warnings
from pathlib import Path


def extract_FEFF(path):
    """Extract spectrum from FEFF output directory.
    Parameters
    ----------
    path : str
        FEFF output directory.
    Returns
    -------
    dict_output
        A dictionary containing the spectrum.
        ``{"label": "...", "energy": [...], "spectrum": [...], "efermi": ...,
        "normalization_constant": ...}``
    """
    dict_output = {}
    dict_output["label"] = "cross-section"
    # FEFF outputs the X-ray absorption cross-section.
    # This label may be used for later analysis.

    path = Path(path)
    with open(path / "xmu.dat", "r") as f:
        spectrum = np.loadtxt(f.readlines())
        dict_output["energy"] = spectrum[:, 0]
        dict_output["spectrum"] = spectrum[:, 3]
        f.seek(0)
        dict_output["normalization_constant"] = [
            float(line.split()[-1])
            for line in f.readlines()
            if "normalize mu" in line
        ][0]
    # Extract sepctrum and normalization constant
    # The normalization constant can be used to normalize feff to
    # absolute cross section in the unit of Angstrom^2

    efermi = spectrum[spectrum[:, 2] == 0, 1]
    # Extract Fermi energy at k=0

    dict_output["efermi"] = efermi[0]

    return dict_output


def extract_VASP(path):
    """Extract unbroadened raw spectrum from VASP output directory.
    Parameters
    ----------
    path : str
        VASP output directory containing ``OUTCAR`` file.

    Returns
    -------
    dict_output
        A dictionary containing the spectrum. The spectrum is averaged from
        all polarization directions.
        The volume of the super cell is extracted in the unit of Angstrom^3.
        Fermi energy and total energy are also extracted if these are contained
        in OUTCAR, in the unit of eV.
        ``{"label": ..., "energy": [...], "efermi": ..., "xx", [...],
        "yy", [...], "zz", [...], "spectrum": [...], "volume": ...,
        "total_energy": ...}``
    """
    dict_output = {}
    dict_output["label"] = "epsilon2"  # VASP outputs the imaginary part of
    # the dielectric tensor. This label may be used for later analysis.

    path = Path(path)
    file_outcar = path / "OUTCAR"
    volume = None
    with open(file_outcar, "r") as f:
        for idx, line in enumerate(f.readlines()):
            if "volume of cell" in line and volume is None:
                volume = float(line.split()[4])
            if "Fermi energy:" in line:
                efermi = float(line.split()[2])
            if "TOTEN" in line:
                total_energy = float(line.split()[-2])
            if "IMAGINARY DIELECTRIC FUNCTION" in line:
                idx_start = idx
            if "REAL DIELECTRIC FUNCTION" in line:
                idx_end = idx
                break
        f.seek(0)
        lines_dielectric = f.readlines()[(idx_start + 3) : (idx_end - 1)]
    mu = np.array([xx.split() for xx in lines_dielectric]).astype(float)
    dict_output["energy"] = mu[:, 0]
    dict_output["efermi"] = efermi
    dict_output["xx"] = mu[:, 1]
    dict_output["yy"] = mu[:, 2]
    dict_output["zz"] = mu[:, 3]
    # xx, yy, zz are the three polarization directions

    dict_output["spectrum"] = mu[:, 1:4].mean(axis=1)
    dict_output["volume"] = volume
    dict_output["total_energy"] = total_energy

    return dict_output


def extract_XSpectra(path, es_out_file=None):
    """Extract spectrum from XSpectra output directory.
    Parameters
    ----------
    path : str
        XSpectra output directory.

    Returns
    -------
    dict_output
        A dictionary containing the spectrum. The spectrum is averaged from
        all polarization directions in ``path``
        Fermi level (in eV) and Total energy (in eV) are extracted if the
        path and file name of es.out file is given.
        ``{"label": ..., "energy": [...], "spectrum": [...], "dipole1", [...],
        "dipole2", [...], "dipole3", [...], "efermi": ...,
        "total_energy": ...}``

    Raises
    ------
    AssertionError
        If the output path does not contain ``xanes.dat`` file.
    """
    dict_output = {}
    dict_output["label"] = "cross-section"
    # XSpectra outputs the X-ray absorption cross-section.
    # This label may be used for later analysis.

    path = Path(path)
    num_polar = 0
    efermi = None
    spectra = 0
    for num_polar, sub_path in enumerate(path.rglob("xanes.dat"), start=1):
        if efermi is None:
            file_out = Path("/".join(Path(sub_path).parts[:-1])) / "xanes.out"
            if file_out.exists():
                with open(file_out, "r") as f:
                    dict_output["efermi"] = [
                        float(line.split()[2])
                        for line in f.readlines()
                        if "ef    [eV]" in line
                    ][0]
        # Try to extract Fermi energy from ``xanes.out`` file.

        if "energy" not in dict_output:
            energy_spectrum = np.loadtxt(sub_path, usecols=(0, 1))
            dict_output["energy"] = energy_spectrum[:, 0]
            spectrum = energy_spectrum[:, 1]
        else:
            spectrum = np.loadtxt(sub_path, usecols=(1))
        dict_output[sub_path.parts[-2]] = spectrum
        spectra += spectrum
    assert num_polar > 0, "'xanes.dat' file not found in %s" % path
    # Raise error if there is no ``xanes.dat`` file in the directory.

    if num_polar > 3:
        warnings.warn(
            "More than three 'xanes.dat' files in %s. \
        Might be a problem."
            % path
        )
    # Print a reminder if more than three polarizations in the directory.

    dict_output["spectrum"] = spectra / num_polar
    if es_out_file is not None:
        with open(es_out_file, "r") as f:
            dict_output["total_energy"] = [
                float(line.split()[-2])
                * physical_constants["Rydberg constant times hc in eV"][0]
                for line in f.readlines()
                if "!" in line
            ][0]
    # Extract total energy if filename of es.out file is given

    return dict_output


def extract_OCEAN(path, scf_out_file=None):
    """Extract spectrum from OCEAN output directory.
    Parameters
    ----------
    path : str
        OCEAN spectra output directory.
    scf_out_file: str, optional
        Path of the OCEAN output file ``scf.out``.

    Returns
    -------
    dict_output
        A dictionary containing the spectra from different sites.
        Each spectrum is averaged from all polarization directions in ``path``.
        Fermi level and total energy are extracted if ``scf.out`` file is
        provided, with unit eV
        ``{
            "Ti": {
                "0001_1s": {"label": ..., "energy": [...], "spectrum": [...],
                "01": [...], "02": [...], "03": [...], "efermi": ...,
                "total_energy": ...},
                "0002_1s": {...},
                ...
                },
           }``
        The sub-dictionary of one site is dict_output["Ti"]["0001_1s"], which
        has the same structure and keys as the output of the parsing functions
        for other codes.

    Raises
    ------
    AssertionError
        If the output path does not contain any file that starts with
        ``absspct``.
    """
    dict_output = {}

    efermi = None
    total_energy = None
    if scf_out_file is not None:
        with open(scf_out_file, "r") as f:
            for line in f.readlines():
                if "Fermi energy" in line:
                    efermi = float(line.split()[4])
                if "!" in line:
                    total_energy = (
                        float(line.split()[-2])
                        * physical_constants["Rydberg constant times hc in eV"][
                            0
                        ]
                    )
    # Try to extract Fermi energy from ``scf.out`` file.

    path = Path(path)
    num_spectra = 0
    for num_spectra, sub_path in enumerate(path.rglob("absspct*"), start=1):
        element = sub_path.parts[-1].split(".")[0].split("_")[1]
        absorber = (
            sub_path.parts[-1].split(".")[1].split("_")[0]
            + "_"
            + sub_path.parts[-1].split(".")[1].split("_")[1]
        )
        polarization = sub_path.parts[-1].split(".")[1].split("_")[-1]
        if element not in dict_output:
            dict_output[element] = {}
        if absorber not in dict_output[element]:
            energy_spectrum = np.loadtxt(sub_path, usecols=(0, 2))
            spectrum = energy_spectrum[:, 1]
            dict_output[element][absorber] = {
                "label": "epsilon2",
                # OCEAN outputs the imaginary part of the dielectric tensor.
                # This label may be used for later analysis.
                "energy": energy_spectrum[:, 0],
                "spectrum": 0,
                # The spectrum value will be assigned later
                "num_polar": 1,
            }
            if efermi is not None:
                dict_output[element][absorber]["efermi"] = efermi
            if total_energy is not None:
                dict_output[element][absorber]["total_energy"] = total_energy
        else:
            spectrum = np.loadtxt(sub_path, usecols=(2))
            dict_output[element][absorber]["num_polar"] += 1
        dict_output[element][absorber][polarization] = spectrum
        dict_output[element][absorber]["spectrum"] += spectrum
    assert num_spectra > 0, "'absspct...' file not found in %s" % path
    # Raise error if there is no ``absspct...`` file in the directory.

    for element in dict_output:
        for absorber in dict_output[element]:
            dict_output[element][absorber]["spectrum"] /= dict_output[element][
                absorber
            ]["num_polar"]
            del dict_output[element][absorber]["num_polar"]
    return dict_output


def extract_exciting(path, INFO_out_file=None):
    """Extract spectrum from exciting output directory.

    Parameters
    ----------
    path : str
        exciting spectra output directory.
    INFO_out_file: str, optional
        Path of the exciting output file ``INFO.out``.

    Returns
    -------
    dict_output
        A dictionary containing the spectrum. The spectrum is averaged from all
        polarization directions in ``path``.
        Fermi level and total energy are extracted if ``INFO.OUT`` file is
        provided, with unit eV.
        ``{"label": ..., "energy": [...], "spectrum": [...], "11": [...],
        "22": [...], "33": [...], "efermi": ..., "total_energy": ...}``

    Raises
    ------
    AssertionError
        If the output path does not contain any file that starts with
        ``EPSILON`` and ends with ``OUT``.
    """

    dict_output = {}
    dict_output["label"] = "epsilon2"
    # EXCITING outputs the imaginary part of the dielectric tensor.
    # This label may be used for later analysis.

    path = Path(path)

    num_polar = 0
    spectra = 0
    for num_polar, sub_path in enumerate(path.rglob("EPSILON*OUT"), start=1):
        if "energy" not in dict_output:
            energy_spectrum = np.loadtxt(sub_path, usecols=(0, 2))
            dict_output["energy"] = energy_spectrum[:, 0]
            spectrum = energy_spectrum[:, 1]
        else:
            spectrum = np.loadtxt(sub_path, usecols=(2))
        polarization = "".join(
            [
                ii
                for ii in sub_path.parts[-1].split(".")[0].split("_")[-1]
                if ii.isdigit()
            ]
        )
        dict_output[polarization] = spectrum
        spectra += spectrum
    assert num_polar > 0, "'EPSILON...OUT' file not found in %s" % path
    # Raise error if there is no ``EPSILON...OUT`` file in the directory.

    if num_polar > 3:
        warnings.warn(
            "More than three 'EPSILON...OUT' files in %s. \
        Might be a problem."
            % path
        )
    # Print a reminder if more than three polarizations in the directory.

    dict_output["spectrum"] = spectra / num_polar

    if INFO_out_file is not None:
        lines_fermi = []
        with open(INFO_out_file, "r") as f:
            for line in f.readlines():
                if "Fermi energy" in line:
                    lines_fermi.append(line)
                if "Total energy" in line:
                    total_energy = (
                        float(line.split()[3])
                        * physical_constants["Hartree energy in eV"][0]
                    )
            dict_output["efermi"] = (
                float(lines_fermi[-2].split()[3])
                * physical_constants["Hartree energy in eV"][0]
            )
            dict_output["total_energy"] = total_energy
    # Extract Fermi energy and total energy from ``INFO.OUT`` file.

    return dict_output

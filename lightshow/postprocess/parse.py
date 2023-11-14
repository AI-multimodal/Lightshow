import numpy as np
import os
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
        The Fermi level is extracted if ``feff.out`` file is in the directory.
        ``{"spectrum": [..., ...], "efermi": ...}``

    Raises
    ------
    OSError
        If the output path does not contain ``xmu.dat`` file. 
    """
    
    dict_output = {}
    
    path = Path(path)
    try: 
        spectrum = np.loadtxt(path / "xmu.dat")
    except OSError as err: 
        print(err)
        return
    dict_output['spectrum'] = np.column_stack((spectrum[:,0],spectrum[:,3]))
    
    ### Try to extract Fermi energy from ``feff.out`` file
    file_out = path / "feff.out"
    lines_fermi = []
    try: 
        with open(file_out, 'r') as f:
            for line in f.readlines():
                if 'New Fermi level' in line:
                    lines_fermi.append(line)
        try: 
            dict_output['efermi'] = float(lines_fermi[-1].split()[4])
        except NameError as err: 
            print('Cannot find Fermi energy in', file_out)                    
    except: 
        print(file_out, 'file not found. Fermi energy not extracted.')

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
        A dictionary containing the spectrum. The spectrum is averaged from all polarization directions. 
        The volume of the super cell and fermi energy are also extracted if these are contained in OUTCAR. 
        ``{"spectrum": [..., ...], "volume": ..., "efermi": ...}``

    Raises
    ------
    OSError
        If the output path does not contain ``OUTCAR`` file. 
    """
    
    dict_output = {}
    
    path = Path(path)
    file_outcar = path / "OUTCAR"
    volume=None
    try: 
        with open(file_outcar, 'r') as f:
            for idx, line in enumerate(f.readlines()):
                if 'volume of cell' in line and volume is None: 
                    volume = float(line.split()[4])
                if 'Fermi energy:' in line:
                    efermi = float(line.split()[2])
                if 'IMAGINARY DIELECTRIC FUNCTION' in line:
                    idx_start = idx
                if 'REAL DIELECTRIC FUNCTION' in line: 
                    idx_end = idx
                    break
            f.seek(0)
            try: 
                lines_dielectric = f.readlines()[(idx_start+3):(idx_end-1)]
            except: 
                print('Wrong OUTCAR file! ')
                return
    except OSError as err: 
        print(err)
        return
    mu = np.array([xx.split() for xx in lines_dielectric]).astype(float)
    mu = np.array([mu[:, 0], mu[:, 1:4].mean(axis=1)]).T
    dict_output['spectrum'] = mu
    dict_output['volume'] = volume
    try: dict_output['efermi'] = efermi 
    except: pass

    return dict_output


def extract_XSpectra(path): 
    """Extract spectrum from XSpectra output directory. 
    
    Parameters
    ----------
    path : str
        XSpectra output directory. 

    Returns
    -------
    dict_output
        A dictionary containing the spectrum. The spectrum is averaged from all polarization directions in ``path``
        Fermi level is extracted if ``xanes.out`` file is in the directory.
        ``{"spectrum": [..., ...], "efermi": ...}``

    Raises
    ------
    AssertionError
        If the output path does not contain ``xanes.dat`` file. 
    """
    
    dict_output = {}
    
    path = Path(path)
    num_polar = 0
    spectra = None
    efermi = None
    for sub_path in path.rglob('xanes.dat'): 
        num_polar += 1
        if spectra is None: 
            spectra = np.loadtxt(sub_path, usecols=(0,1))
        else: 
            spectra[:,1] += np.loadtxt( sub_path, usecols=(1)  )
        ### Try to extract Fermi energy from ``xanes.out`` file. 
        if efermi is None: 
            file_out = Path("/".join(Path(sub_path).parts[:-1]))/"xanes.out"
            try:
                with open(file_out, 'r') as f:
                    for line in f.readlines():
                        if 'ef    [eV]' in line:
                            efermi = float(line.split()[2])
                            dict_output['efermi'] = efermi
            except: 
                print(file_out, 'not found. Fermi energy not extracted.')
    
    ### Raise error if there is no ``xanes.dat`` file in the directory.                         
    try: 
        assert num_polar > 0, "\'xanes.dat\' file not found in %s"%path 
    except AssertionError as msg: 
        print(msg)
        return
    ### Print a reminder if there are more than three polarizations in the directory. 
    try: 
        assert num_polar <= 3, "More than three \'xanes.dat\' files in %s. Might be a problem."%path
    except AssertionError as msg: 
        print(msg)
        
    spectra[:,1] /= num_polar
    dict_output['spectrum'] = spectra

    return dict_output


def extract_OCEAN(path, scf_path=None): 
    """Extract spectrum from OCEAN output directory. 
    
    Parameters
    ----------
    path : str
        OCEAN spectra output directory. 
    scf_path: str, optional
        OCEAN output directory where ``scf.out`` file is stored. 
        If scf_path is not specified, will use the spectra path. 

    Returns
    -------
    dict_output
        A dictionary containing the spectra from different sites. Each spectrum is averaged from all polarization directions in ``path``.
        Fermi level is extracted if ``scf.out`` file is in the ``scf_path`` directory.
        ``{
            "Ti": {"0001_1s": {"spectrum": [..., ...] }, "0002_1s": {"spectrum": [..., ...] } }, 
            "efermi": ...
           }``

    Raises
    ------
    AssertionError
        If the output path does not contain any file that starts with ``absspct``. 
    """
    
    if scf_path is None:
        scf_path = path
    
    dict_output = {}
    
    path = Path(path)
    num_spectra = 0
    spectra = None
    for sub_path in path.rglob('absspct*'): 
        num_spectra += 1
        element = sub_path.parts[-1].split('.')[0].split('_')[1]
        absorber = sub_path.parts[-1].split('.')[1].split('_')[0] + '_' + sub_path.parts[-1].split('.')[1].split('_')[1]
        if element not in dict_output.keys(): 
            dict_output[element] = {}
        else: 
            if absorber not in dict_output[element].keys(): 
                dict_output[element][absorber] = {
                    'spectrum': np.loadtxt( sub_path, usecols=(0,2)  ), 
                    'num_polar': 1
                }
            else: 
                dict_output[element][absorber]['spectrum'][:,1] += np.loadtxt( sub_path, usecols=(2)  )
                dict_output[element][absorber]['num_polar'] += 1

    ### Raise error if there is no ``absspct...`` file in the directory.                     
    try: 
        assert num_spectra > 0, "\'absspct...\' file not found in %s"%path 
    except AssertionError as msg: 
        print(msg)
        return
    
    for element in dict_output.keys(): 
        for absorber in dict_output[element].keys():
            dict_output[element][absorber]['spectrum'][:,1] /= dict_output[element][absorber]['num_polar']
            del dict_output[element][absorber]['num_polar']

    ### Try to extract Fermi energy from ``scf.out`` file. 
    file_out = Path(scf_path) / 'scf.out'
    try: 
        with open(file_out, 'r') as f:
            for line in f.readlines():
                if 'Fermi energy' in line:
                    efermi = float(line.split()[4])
        try: 
            dict_output['efermi'] = efermi
        except NameError as err: 
            print('Cannot find Fermi energy in', file_out)
    except OSError: 
        print(file_out, 'not found. Fermi energy not extracted.')
    
    return dict_output


def extract_exciting(path, INFO_path=None): 
    """Extract spectrum from exciting output directory. 
    
    Parameters
    ----------
    path : str
        exciting spectra output directory. 
    INFO_path: str, optional
        exciting output directory where ``INFO.OUT`` file is stored. 
        If INFO_path is not specified, will use the spectra path. 
        
    Returns
    -------
    dict_output
        A dictionary containing the spectrum. The spectrum is averaged from all polarization directions in ``path``.
        Fermi level is extracted if ``INFO.OUT`` file is in the ``INFO_path`` directory.
        ``{"spectrum": [..., ...], "efermi": ...}``

    Raises
    ------
    AssertionError
        If the output path does not contain any file that starts with ``EPSILON`` and ends with ``OUT``. 
    """
    
    if INFO_path is None:
        INFO_path = path
    dict_output = {}
    
    path = Path(path)
    num_polar = 0
    spectra = None
    for sub_path in path.rglob('EPSILON*OUT'): 
        num_polar += 1
        if spectra is None: 
            spectra = np.loadtxt(sub_path, usecols=(0,2))
        else: 
            spectra[:,1] += np.loadtxt( sub_path, usecols=(2)  )

    ### Raise error if there is no ``EPSILON...OUT`` file in the directory.                     
    try: 
        assert num_polar > 0, "\'EPSILON...OUT\' file not found in %s"%path 
    except AssertionError as msg: 
        print(msg)
        return
    ### Print a reminder if there are more than three polarizations in the directory. 
    try: 
        assert num_polar <= 3, "More than three \'EPSILON...OUT\' files in %s. Might be a problem."%path
    except AssertionError as msg: 
        print(msg)
        
    spectra[:,1] /= num_polar
    dict_output['spectrum'] = spectra

    ### Try to extract Fermi energy from ``INFO.OUT`` file. 
    file_out = Path(INFO_path) / 'INFO.OUT'
    lines_fermi = []
    try: 
        with open(file_out, 'r') as f:
            for line in f.readlines():
                if 'Fermi energy' in line:
                    lines_fermi.append(line)
        try: 
            dict_output['efermi'] = float(lines_fermi[-2].split()[3])
        except: 
            print('Cannot find Fermi energy in', file_out) 
    except: 
        print(file_out, 'not found. Fermi energy not extracted.')

    return dict_output



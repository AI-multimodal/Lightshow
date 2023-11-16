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
        ``{"spectrum": [..., ...], "efermi": ..., "normalization_constant": ...}``

    Raises
    ------
    OSError
        If the output path does not contain ``xmu.dat`` file. 
    """
    
    dict_output = {}    
    path = Path(path)
    
    # Extract sepctrum 
    try: 
        spectrum = np.loadtxt(path / "xmu.dat")
    except OSError as err: 
        print(err)
        return
    dict_output['spectrum'] = np.column_stack((spectrum[:,0],spectrum[:,3]))

    # Extract Fermi energy at k=0
    efermi = spectrum[spectrum[:,2]==0,1]
    try: 
        dict_output['efermi'] = efermi[0]
    except: 
        print("Failed to extract Fermi energy.")
        
    # Extract the normalization constant from the header
    # This number can be used to normalize feff to absolute cross section in the unit of Angstrom^2 
    with open(path/"xmu.dat", 'r') as f:
        for line in f.readlines():
            if 'normalize mu' in line:
                normalization_constant = float(line.split()[-1])
                dict_output['normalization_constant'] = normalization_constant
                break

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
        The volume of the super cell is extracted in the unit of Angstrom^3.
        Fermi energy and total energy are also extracted if these are contained in OUTCAR, in the unit of eV. 
        ``{"spectrum": [..., ...], "volume": ..., "efermi": ..., "total_energy": ...}``

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
                if 'TOTEN' in line:
                    total_energy = float(line.split()[-2])
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
    try: dict_output['total_energy'] = total_energy 
    except: pass

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
        A dictionary containing the spectrum. The spectrum is averaged from all polarization directions in ``path``
        Fermi level (in eV) is extracted if ``xanes.out`` file is in the directory.
        Total energy (in eV) is extracted if the path and file name of es.out file is given. 
        ``{"spectrum": [..., ...], "efermi": ..., "total_energy": ...}``

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
                            break
                try: 
                    dict_output['efermi'] = efermi
                except NameError as err: 
                    print('Cannot find Fermi energy in', file_out)                            
            except: 
                print(file_out, 'not found. Fermi energy not extracted.')
    # Raise error if there is no ``xanes.dat`` file in the directory.                         
    try: 
        assert num_polar > 0, "\'xanes.dat\' file not found in %s"%path 
    except AssertionError as msg: 
        print(msg)
        return
    # Print a reminder if there are more than three polarizations in the directory. 
    try: 
        assert num_polar <= 3, "More than three \'xanes.dat\' files in %s. Might be a problem."%path
    except AssertionError as msg: 
        print(msg)
    spectra[:,1] /= num_polar
    dict_output['spectrum'] = spectra
    
    # Extract total energy if filename of es.out file is given
    if es_out_file is not None: 
        try: 
            with open(es_out_file, 'r') as f:
                for line in f.readlines():
                    if '!' in line:
                        total_energy = float(line.split()[-2]) * 13.605703976
                        break
            try: 
                dict_output['total_energy'] = total_energy
            except NameError as err: 
                print('Cannot find total energy in', es_out_file)                        
        except: 
            print(es_out_file, 'not found. Total energy not extracted.')                

    return dict_output


def extract_OCEAN(path, scf_out_file=None): 
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
        Fermi level and total energy are extracted if ``scf.out`` file is provided, with unit eV
        ``{
            "Ti": {"0001_1s": {"spectrum": [..., ...] }, "0002_1s": {"spectrum": [..., ...] } }, 
            "efermi": ...,
            "total_energy": ...
           }``

    Raises
    ------
    AssertionError
        If the output path does not contain any file that starts with ``absspct``. 
    """
    
    dict_output = {}
    
    path = Path(path)
    if scf_out_file is None:
        scf_out_file = path / "scf.out"
        
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

    # Raise error if there is no ``absspct...`` file in the directory.                     
    try: 
        assert num_spectra > 0, "\'absspct...\' file not found in %s"%path 
    except AssertionError as msg: 
        print(msg)
        return
    
    for element in dict_output.keys(): 
        for absorber in dict_output[element].keys():
            dict_output[element][absorber]['spectrum'][:,1] /= dict_output[element][absorber]['num_polar']
            del dict_output[element][absorber]['num_polar']

    # Try to extract Fermi energy from ``scf.out`` file. 
    file_out = scf_out_file
    try: 
        with open(file_out, 'r') as f:
            for line in f.readlines():
                if 'Fermi energy' in line:
                    efermi = float(line.split()[4])
                if '!' in line:
                    total_energy = float(line.split()[-2]) * 13.605703976
        try: 
            dict_output['efermi'] = efermi
        except NameError as err: 
            print('Cannot find Fermi energy in', file_out)
        try: 
            dict_output['total_energy'] = total_energy
        except: 
            print('Cannot find total energy in', file_out) 
    except OSError: 
        print(file_out, 'not found. Fermi energy and total energy not extracted.')
    
    return dict_output


def extract_exciting(path, INFO_out_file=None): 
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
        Fermi level and total energy are extracted if ``INFO.OUT`` file is provided, with unit eV.
        ``{"spectrum": [..., ...], "efermi": ..., "total_energy": ...}``

    Raises
    ------
    AssertionError
        If the output path does not contain any file that starts with ``EPSILON`` and ends with ``OUT``. 
    """
    
    dict_output = {}
    
    path = Path(path)
    if INFO_out_file is None:
        INFO_out_file = path / "INFO.OUT"
    
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
    # Print a reminder if there are more than three polarizations in the directory. 
    try: 
        assert num_polar <= 3, "More than three \'EPSILON...OUT\' files in %s. Might be a problem."%path
    except AssertionError as msg: 
        print(msg)
        
    spectra[:,1] /= num_polar
    dict_output['spectrum'] = spectra

    # Extract Fermi energy and total energy from ``INFO.OUT`` file. 
    file_out = INFO_out_file
    lines_fermi = []
    try: 
        with open(file_out, 'r') as f:
            for line in f.readlines():
                if 'Fermi energy' in line:
                    lines_fermi.append(line)
                if 'Total energy' in line:
                    total_energy = float(line.split()[3])*27.211407953
        try: 
            dict_output['efermi'] = float(lines_fermi[-2].split()[3])*27.211407953
        except: 
            print('Cannot find Fermi energy in', file_out) 
        try: 
            dict_output['total_energy'] = total_energy
        except: 
            print('Cannot find total energy in', file_out) 
    except: 
        print(file_out, 'not found. Fermi energy and total energy not extracted.')

    return dict_output



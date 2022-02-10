# This is a simple tool to compare to plots
# Fanchen Meng, 2022
# FC TODO 1. paser for different codes
#         2. compatible with convergece

from xanes_bench.CurveComparisons.compare_pasers import XSplot, OCEANplot, Excitingplot
from xanes_bench.CurveComparisons.compare_utils import comparePlots
from pathlib import Path
import numpy as np
from collections import OrderedDict
import json
import matplotlib.pyplot as plt
from pymatgen.core import Structure
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
import sys


def confirm_folders(mpid, code):
    k_folder = Path(mpid) / code / 'k.txt'
    k_list = np.loadtxt(k_folder, usecols=[0,1,2,5])
    spectra_folder_dict = OrderedDict()
    for k in k_list:
        spectra_name = 'Spectra-' + str(int(k[0])) + '-' + str(int(k[1])) + '-' + str(int(k[2]))
        spectra_folder = Path(mpid) / code / spectra_name
        if spectra_folder.is_dir():
            spectra_folder_dict[spectra_name] = k[-1]
    return spectra_folder_dict

def find_sites(mpid, target):
    sites = []
    json_file = json.load(open(f'{mpid}/{mpid}.json', 'r'))
    struct = Structure.from_dict(json_file)
    space_group_analyzer = SpacegroupAnalyzer(struct)
    symm_struct = space_group_analyzer.get_symmetrized_structure()
    
    for idx in symm_struct.equivalent_indices:
        if str(symm_struct[idx[0]].specie) == target:
            sites.append(str(idx[0]))
    return sites


def calculate_similarity_last(mpid, code):
    site_list = find_sites(mpid, target = 'Ti')
    spectra_folder_dict = confirm_folders(mpid, code)
    num_folders = len(spectra_folder_dict)
    
    last_spectra_fn = Path(mpid) / code / list(spectra_folder_dict.keys())[num_folders-1]
    if code == 'XS':
        last_spectra = XSplot(last_spectra_fn, site_list)
    elif code == 'OCEAN':
        last_spectra = OCEANplot(last_spectra_fn, site_list)
    elif code == 'EXCITIG':
        last_spectra = Excitingplot(last_spectra_fn, site_list)

    while not last_spectra.exists():
        num_folders = num_folders - 1
        last_spectra_fn = Path(mpid) / code / list(spectra_folder_dict.keys())[num_folders-1]
        if code == 'XS':
            last_spectra = XSplot(last_spectra_fn, site_list)
        elif code == 'OCEAN':
            last_spectra = OCEANplot(last_spectra_fn, site_list)
        elif code == 'EXCITIG':
            last_spectra = Excitingplot(last_spectra_fn, site_list)
    print(f'INFO: Most Converged Spectra is in {last_spectra_fn}')

    out = dict()
    for i in spectra_folder_dict:
        this_spectra_fn = Path(mpid) / code / i
        if code == 'XS':
            this_spectra = XSplot(this_spectra_fn, site_list)
        elif code == 'OCEAN':
            this_spectra = OCEANplot(this_spectra_fn, site_list)
        elif code == 'EXCITIG':
            this_spectra = Excitingplot(this_spectra_fn, site_list)
        if this_spectra.exists():
            spearman = comparePlots(0, this_spectra.spectra, last_spectra.spectra)[1]
            out[spectra_folder_dict[i]] = np.log10(1-spearman)
        else:
            continue
    return out

def calculate_similarity_next(mpid, code):
    site_list = find_sites(mpid, target = 'Ti')
    spectra_folder_dict = confirm_folders(mpid, code)
    
    out = dict()
    spectra_folder_list = list(spectra_folder_dict.keys()) 
    
    for i in range(len(spectra_folder_list)-1):
        
        first_spectra_fn = Path(mpid) / code / spectra_folder_list[i]
        second_spectra_fn = Path(mpid) / code / spectra_folder_list[i+1]
        
        if code == 'XS':
            first_spectra = XSplot(first_spectra_fn, site_list)
            second_spectra = XSplot(second_spectra_fn, site_list)
        elif code == 'OCEAN':
            first_spectra = OCEANplot(first_spectra_fn, site_list)
            second_spectra = OCEANplot(second_spectra_fn, site_list)
        elif code == 'EXCITING':
            first_spectra = Excitingplot(first_spectra_fn, site_list)
            second_spectra = Excitingplot(second_spectra_fn, site_list)
            
        spearman = comparePlots(0, first_spectra.spectra, second_spectra.spectra)[1]
        out[spectra_folder_dict[spectra_folder_list[i]]] = np.log10(1-spearman)
    return out

def main():
    
    
    if sys.argv[1] == 'x':
        code = 'XS'
    elif sys.argv[1] == 'o':
        code = 'OCEAN'
    elif sys.argv[1] == 'e':
        code = 'EXCITING'

    mpid = sys.argv[2]

    out_last = calculate_similarity_last(mpid, code)
    out_next = calculate_similarity_next(mpid, code)
    with open(f'spearman_last_{code}_{mpid}.dat', 'w') as f:
        print(f'INFO: comparing with next spectra')
        print(f'INFO: effective radius (Bohr)    log(1-Spearman)')
        for i,j in out_last.items():
            print(i, j)
            txt = str(i) + '  ' + str(j) + '\n'
            f.write(txt)

    with open(f'spearman_next_{code}_{mpid}.dat', 'w') as f:
        print(f'INFO: comparing with next spectra')
        print(f'INFO: effective radius (Bohr)    log(1-Spearman)')
        for i,j in out_next.items():
            print(i, j)
            txt = str(i) + '  ' + str(j) + '\n'
            f.write(txt)


if __name__ == '__main__':
    main()


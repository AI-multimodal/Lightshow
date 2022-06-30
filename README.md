# xanes_bench
- Benchmark code and data for multiple simulation method
- Create input files for qe/xspectra and ocean based on a material ID from the Materials Project

## Quickstart
1. Put your materials project API key in the file mp.key ( https://materialsproject.org/dashboard )
2. Run "python fetchSingle.py mp-390" (will run anatase)


## Limitations/To-do list: (not exhaustive)
1. Fix up symmetry dipole/quad options
2. JSON files with convergence options for xspectra/ocean (remove from script)
3. Unify core-hole broadening (ocean uses the FEFF look-up table)

# Installation Prerequisites
* For XSpectra, need to use latest version of pymatgen. [See here for more info](https://github.com/materialsproject/pymatgen/pull/2178). 

# Install
```console
conda create --name xas_ben python=3
conda activate xas_ben
git clone git@github.com:xhqu1981/xanes_bench.git
cd xanes_bench
pip install -r requirements.txt
python setup.py develop
```
# Test systems
| MP ID | Name |
|---------|---------|
|mp-2657 | rutile TiO2|
|mp-390 | anatase TiO2|
|mp-1840 | brookite TiO2|
|mp-1203 | TiO |

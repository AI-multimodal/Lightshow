# xanes_bench
- Benchmark code and data for multiple simulation method
- Create input files for qe/xspectra and ocean based on a material ID from the Materials Project

## Quickstart
1. Put your materials project API key in the file mp.key ( https://materialsproject.org/dashboard )
2. Overwrite the defaults ase/io/espresso.py file with the one in Xspectra
  - Either: copy it into the installed location of ase ( site-packages/ase/io/epresso.py )
  - Or: copy all of ase locally, and write over ./ase/io/espresso.py
3. Run "python fetchSingle.py 390" (will run anatase)


## Limitations/To-do list: (not exhaustive)
1. Fix up symmetry dipole/quad options
2. Separate xspectra to its own python file
3. JSON files with convergence options for xspectra/ocean (remove from script)
4. Unify core-hole broadening (ocean uses the FEFF look-up table)


# Install
```console
conda create --name xas_ben python=3
conda activate xas_ben
git clone git@github.com:xhqu1981/xanes_bench.git
cd xanes_bench
python setup.py develop
python apply_patch.py --force
```
# Test systems
| MP ID | Name |
|---------|---------|
|mp-2657 | rutile TiO2|
|mp-390 | anatase TiO2|
|mp-1840 | brookite TiO2|
|mp-1203 | TiO |

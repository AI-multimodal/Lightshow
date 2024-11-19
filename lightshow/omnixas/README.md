# Omnixas Integration in Lightshow

## Description

Self-contained modules to predcit site-XAS [`pymatgen`](https://pymatgen.org/) structures using Omnixas models.

## Models

- `v1.1.1`: Available in [Omnixas Repository](https://github.com/AI-multimodal/OmniXAS/tree/main/models)

## Input

- Input is `pymatgen` `Structure` object. Useful I/O methods for `pymatgen` structures are documented in [pyamtgen documentation](https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.IStructure) which includes `from_file`, `from_id`, `from_sites`, `from_spacegroup`, `from_str`.

## Supported File Formats

- use [`from_file`](https://pymatgen.org/pymatgen.core.html#pymatgen.core.structure.IMolecule.from_file) to load `CIF`, `POSCAR`, `CONTCAR`, `CHGCAR`, `LOCPOT`, `vasprun.xml`, `CSSR`, `Netcdf` and `pymatgen's JSON-serialized structures`.

## Usage

```python
import matplotlib.pyplot as plt
material_structure_file = "mp-1005792/POSCAR"
strucutre = PymatgenStructure.from_file(material_structure_file)
spectrum = XASModel(element="Cu", type="FEFF").predict(strucutre, 8)
plt.plot(spectrum)
```

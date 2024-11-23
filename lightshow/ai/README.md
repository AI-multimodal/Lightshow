# X-ray Absorption Spectroscopy Prediction Integration in Lightshow

> [!WARNING]
> This feature is a work in progress! Expect breaking changes!

Using [machine-learned neural networks](https://doi.org/10.48550/arXiv.2409.19552) stacked on top of the [M3GNet architecture](https://github.com/materialsvirtuallab/m3gnet), we have built in direct XANES spectrum prediction. To use this feature, you should install Lightshow's AI dependencies:
```bash
pip install "lightshow[ai]"
```
and then use

```python
from pymatgen.core.structure import Structure
from lightshow.ai.models import predict

structure = Structure.from_file(...)  # or however
predictions = predict(structure, "Ti", "FEFF")
```

where `predictions` is a dictionary of the spectra indexed by absorption site.


## Lightshow server

To support the possibility of loading models only once, or putting the machine learning component of this functionality on a separate machine, Lightshow supports running a server separately.

On the server machine:

```bash
pip install "lightshow[ai,server]"
lightshow-serve
```
On the "local" machine:

```bash
python lightshow/ai/test_server.py PATH_TO_STRUCTURE
```

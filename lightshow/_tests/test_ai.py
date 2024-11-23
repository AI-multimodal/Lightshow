from pathlib import Path

from pymatgen.core.structure import Structure

from lightshow.ai.models import XASModel

# Helper testing files
ai_structure_path = (
    Path(__file__).parent.resolve()
    / "ai_examples"
    / "material"
    / "mp-1005792"
    / "POSCAR"
)


def test_predict_api():
    struct = Structure.from_file(ai_structure_path)
    XASModel(element="Cu", spectroscopy_type="FEFF").predict(struct)

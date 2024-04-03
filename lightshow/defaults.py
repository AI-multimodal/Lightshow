"""Helper module for making the API a bit cleaner when dealing with default
parameters."""

from lightshow.parameters.exciting import EXCITING_DEFAULT_CARDS
from lightshow.parameters.feff import FEFF_DEFAULT_CARDS
from lightshow.parameters.ocean import OCEAN_DEFAULT_CARDS
from lightshow.parameters.vasp import (
    VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL,
    VASP_INCAR_DEFAULT_NEUTRAL_POTENTIAL,
    VASP_POTCAR_DEFAULT_ELEMENT_MAPPING,
)
from lightshow.parameters.xspectra import XSPECTRA_DEFAULT_CARDS

__all__ = [
    "FEFF_DEFAULT_CARDS",
    "VASP_INCAR_DEFAULT_NEUTRAL_POTENTIAL",
    "VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL",
    "VASP_POTCAR_DEFAULT_ELEMENT_MAPPING",
    "OCEAN_DEFAULT_CARDS",
    "XSPECTRA_DEFAULT_CARDS",
    "EXCITING_DEFAULT_CARDS",
]

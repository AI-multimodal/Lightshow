"""Helper module for making the API a bit cleaner when dealing with default
parameters."""

# Note that pre-commit + .flake8 ignores this file since it's just a bunch of
# imports...

from lightshow.parameters.vasp import (
    VASP_INCAR_DEFAULT_NEUTRAL_POTENTIAL,
    VASP_INCAR_DEFAULT_COREHOLE_POTENTIAL,
    VASP_POTCAR_DEFAULT_ELEMENT_MAPPING,
)

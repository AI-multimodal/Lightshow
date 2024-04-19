import json
from pathlib import Path

from lightshow.parameters.feff import FEFFParameters


def test_multiplicity_writing(database_from_file, tmp_path):
    target = Path(tmp_path) / "t"
    target.mkdir(exist_ok=True, parents=True)
    R = 10.0
    feff_parameters = FEFFParameters(
        cards={
            "S02": "0",
            "COREHOLE": "RPA",
            "CONTROL": "1 1 1 1 1 1",
            "XANES": "4 0.04 0.1",
            "SCF": "7.0 0 100 0.2 3",
            "FMS": "9.0 0",
            "EXCHANGE": "0 0.0 0.0 2",
            "RPATH": "-1",
        },
        edge="K",
        radius=R - 1.0,
        spectrum="XANES",
        name="FEFF",
    )
    database_from_file.write(
        target,
        absorbing_atoms="all",
        options=[feff_parameters],
        write_unit_cells=True,
        pbar=False,
    )
    for k, v in database_from_file.metadata.items():
        fname = target / k / Path("multiplicity.json")
        with open(fname) as f:
            mult = json.load(f)
        prim_meta = v["primitive"]
        for i_site, n_multi in zip(
            prim_meta["sites"], prim_meta["multiplicities"]
        ):
            i_site = str(i_site)
            assert mult[i_site] == n_multi

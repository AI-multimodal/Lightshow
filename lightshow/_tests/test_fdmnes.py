from pathlib import Path

from pymatgen.core.structure import Structure

from lightshow.parameters.fdmnes import FDMNESParameters


def test_default_parameters(mp_Structure_mp390: Structure, tmp_path: Path):
    # Test on TiO2 O-K default input

    fdmnes_parameters = FDMNESParameters()
    assert fdmnes_parameters.name == "FDMNES"

    write_results = fdmnes_parameters.write(
        tmp_path, structure=mp_Structure_mp390, Z_absorber=8
    )

    input_file = tmp_path / "O_in.txt"
    assert write_results == {
        "pass": True,
        "errors": {},
        "path": str(input_file),
    }
    assert input_file.exists(), "Input file should exist after write operation"

    with open(input_file, "r") as file:
        contents = file.read()
        assert (
            "Spinorbit" in contents
        ), "Spinorbit should be turned on for O-K edge"
        assert (
            "Full_atom" in contents
        ), "Full atom should be turned on for oxides"
        assert "Crystal" in contents, "Crystal should be written in the file"


def test_customized_parameters(mp_Structure_mp390: Structure, tmp_path: Path):
    # Test on customized input on TiO2 Ti-L

    FDMNES_DEFAULT_CARDS = {
        "Energpho": True,
        "Memory_save": True,
        "Quadrupole": False,
        "Relativism": False,
        "Spinorbit": None,
        "SCF": True,
        "SCFexc": False,
        "Screening": False,
        "Full_atom": False,
        "TDDFT": False,
        "PBE96": False,
    }

    cards = FDMNES_DEFAULT_CARDS
    cards["PBE96"] = True

    fdmnes_parameters = FDMNESParameters(cards=cards, edge="L")
    write_results = fdmnes_parameters.write(
        tmp_path, structure=mp_Structure_mp390, Z_absorber=22
    )

    input_file = tmp_path / "Ti_in.txt"
    assert write_results == {
        "pass": True,
        "errors": {},
        "path": str(input_file),
    }
    assert input_file.exists(), "Input file should exist after write operation"

    assert (
        fdmnes_parameters._edge == "L23"
    ), "L edge should be adjust automatically"

    with open(input_file, "r") as file:
        contents = file.read()
        assert "TDDFT" in contents, "TDDFT should be turned on for Ti-L edges"
        assert (
            "Relativism" not in contents
        ), "Relativism should be turned off for this system"

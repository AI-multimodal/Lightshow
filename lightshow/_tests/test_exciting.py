"""
Basic tests for exciting.
"""
import re
from pathlib import Path

from pymatgen.core.structure import Structure
from xml.etree import ElementTree

from lightshow.parameters.exciting import EXCITINGParameters


def test_default_parameters(mp_Structure_mp390: Structure, tmp_path: Path):
    """Test the default parameters of exciting."""
    exciting_parameters = EXCITINGParameters()
    assert exciting_parameters.name == "EXCITING"
    write_results = exciting_parameters.write(
        tmp_path, structure_uc=mp_Structure_mp390, sites=[1]
    )
    assert write_results == {"pass": True, "errors": {}}
    ti_path = tmp_path / "001_Ti"
    assert ti_path.exists(), "Directory for absorbing site was not generated."
    input_file = ti_path / "input.xml"
    assert input_file.exists(), "Input file not found."

    root = ElementTree.parse(input_file).getroot()

    assert root.tag == "input"

    subelements = list(root)
    assert len(subelements) == 4

    title_xml = subelements[0]
    assert title_xml.tag == "title"
    assert title_xml.keys() == []
    assert title_xml.text == "Ti2 O4"

    structure_xml = subelements[1]
    assert structure_xml.tag == "structure"
    assert set(structure_xml.keys()) == {"speciespath", "autormt"}
    assert structure_xml.get("speciespath") == "./", "species path set to ./"
    assert structure_xml.get("autormt") == "true"

    assert (
        len(list(structure_xml)) == 3
    ), "Expect structure tree to have 3 sub-elements"
    crystal_xml = structure_xml[0]
    assert crystal_xml.tag == "crystal", "First subtree is crystal"
    assert crystal_xml.items() == [("scale", "1.8897543760313331")]

    lattice_vectors = list(crystal_xml)
    assert len(lattice_vectors) == 3, "Always expect three lattice vectors"
    assert lattice_vectors[0].items() == [], "Lattice vectors have no items"
    assert (
        lattice_vectors[0].text
        == "     -1.90135500       1.90135500       4.87387600"
    ), "Lattice vector `a` differs from input"
    assert (
        lattice_vectors[1].text
        == "      1.90135500      -1.90135500       4.87387600"
    ), "Lattice vector `b` differs from input"
    assert (
        lattice_vectors[2].text
        == "      1.90135500       1.90135500      -4.87387600"
    ), "Lattice vector `c` differs from input"

    species_ti_xml = structure_xml[1]
    assert species_ti_xml.tag == "species", "Second subtree is species"

    species_o_xml = structure_xml[2]
    assert species_o_xml.tag == "species", "Third subtree is species"

    assert species_ti_xml.items() == [
        ("speciesfile", "Ti.xml")
    ], "species is inconsistent"
    assert species_o_xml.items() == [
        ("speciesfile", "O.xml")
    ], "species is inconsistent"

    atoms_ti = list(species_ti_xml)
    assert len(atoms_ti) == 2, "Number of Ti atoms is wrong"
    assert atoms_ti[0].items() == [
        ("coord", "      0.75000000       0.25000000       0.50000000")
    ], "Coordinate of first Ti differs to input"
    assert atoms_ti[1].items() == [
        ("coord", "      0.00000000       0.00000000       0.00000000")
    ], "Coordinate of second Ti differs to input"

    atoms_o = list(species_o_xml)
    assert len(atoms_o) == 4, "Number of O atoms is wrong"
    assert atoms_o[0].items() == [
        ("coord", "      0.95616300       0.45616300       0.50000000")
    ], "Coordinate of first O differs to input"
    assert re.match(
        r" {6}0\.79383700 {7}0\.79383700 {7}[01]\.00000000",
        atoms_o[1].items()[0][1],
    ), "Coordinate of second O differs to input"
    assert atoms_o[2].items() == [
        ("coord", "      0.54383700       0.04383700       0.50000000")
    ], "Coordinate of third O differs to input"
    assert atoms_o[3].items() == [
        ("coord", "      0.20616300       0.20616300       0.00000000")
    ], "Coordinate of fourth O differs to input"

    groundstate_xml = subelements[2]
    assert groundstate_xml.tag == "groundstate"
    assert (
        len(groundstate_xml.keys()) == 9
    ), "Excpect groundstate to have 9 attributes."
    assert set(groundstate_xml.keys()) == {
        "rgkmax",
        "do",
        "ngridk",
        "xctype",
        "nempty",
        "gmaxvr",
        "lmaxmat",
        "lmaxvr",
        "lmaxapw",
    }, "Excpected groundstate keys differ from found keys."
    assert groundstate_xml.get("rgkmax") == "9.0"
    assert groundstate_xml.get("do") == "fromscratch"
    assert groundstate_xml.get("ngridk") == "3 3 4"
    assert groundstate_xml.get("xctype") == "GGA_PBE"
    assert groundstate_xml.get("nempty") == "200"
    assert groundstate_xml.get("gmaxvr") == "25"
    assert groundstate_xml.get("lmaxmat") == "10"
    assert groundstate_xml.get("lmaxvr") == "10"
    assert groundstate_xml.get("lmaxapw") == "10"

    xs_xml = subelements[3]
    assert xs_xml.tag == "xs"
    assert len(xs_xml.keys()) == 9, "Excpect xs to have 9 attributes."
    assert set(xs_xml.keys()) == {
        "broad",
        "ngridk",
        "xstype",
        "nempty",
        "gqmax",
        "vkloff",
        "tevout",
        "tappinfo",
        "ngridq",
    }, "Excpected xs keys differ from found keys."
    assert xs_xml.get("broad") == "0.0327069"
    assert xs_xml.get("ngridk") == "3 3 4"
    assert xs_xml.get("xstype") == "BSE"
    assert xs_xml.get("nempty") == "150"
    assert xs_xml.get("gqmax") == "4.0"
    assert xs_xml.get("vkloff") == "0.05 0.03 0.13"
    assert xs_xml.get("tevout") == "true"
    assert xs_xml.get("tappinfo") == "true"
    assert xs_xml.get("ngridq") == "3 3 4"

    xs_subelements = list(xs_xml)
    assert len(xs_subelements) == 5, "Expect xs tree to have 5 sub-elements"

    energywindow_xml = xs_subelements[0]
    assert (
        energywindow_xml.tag == "energywindow"
    ), "First xs subelement should be energywindow"
    assert set(energywindow_xml.keys()) == {"intv", "points"}
    assert energywindow_xml.get("intv") == "178.2 180.5"
    assert energywindow_xml.get("points") == "1000"

    screening_xml = xs_subelements[1]
    assert (
        screening_xml.tag == "screening"
    ), "Second xs subelement should be screening"
    assert set(screening_xml.keys()) == {"screentype", "nempty"}
    assert screening_xml.get("screentype") == "full"
    assert screening_xml.get("nempty") == "150"

    bse_xml = xs_subelements[2]
    assert bse_xml.tag == "BSE", "Third xs subelement should be BSE"
    assert set(bse_xml.keys()) == {
        "bsetype",
        "xas",
        "xasedge",
        "nstlxas",
        "distribute",
        "eecs",
        "xasspecies",
        "xasatom",
    }
    assert bse_xml.get("bsetype") == "singlet"
    assert bse_xml.get("xas") == "true"
    assert bse_xml.get("xasedge") == "K"
    assert bse_xml.get("nstlxas") == "1 59"
    assert bse_xml.get("distribute") == "true"
    assert bse_xml.get("eecs") == "1000"
    assert bse_xml.get("xasspecies") == "1"
    assert bse_xml.get("xasatom") == "2"

    plan_xml = xs_subelements[3]
    assert plan_xml.tag == "plan", "Fourth xs subelement should be plan"
    assert plan_xml.items() == []
    doonlys = list(plan_xml)
    assert len(doonlys) == 8, "Excpect 8 plan elements"
    assert doonlys[0].tag == "doonly"
    assert doonlys[0].items() == [("task", "xsgeneigvec")]
    assert doonlys[1].tag == "doonly"
    assert doonlys[1].items() == [("task", "writepmatxs")]
    assert doonlys[2].tag == "doonly"
    assert doonlys[2].items() == [("task", "scrgeneigvec")]
    assert doonlys[3].tag == "doonly"
    assert doonlys[3].items() == [("task", "scrwritepmat")]
    assert doonlys[4].tag == "doonly"
    assert doonlys[4].items() == [("task", "screen")]
    assert doonlys[5].tag == "doonly"
    assert doonlys[5].items() == [("task", "scrcoulint")]
    assert doonlys[6].tag == "doonly"
    assert doonlys[6].items() == [("task", "exccoulint")]
    assert doonlys[7].tag == "doonly"
    assert doonlys[7].items() == [("task", "bse")]

    qpointset_xml = xs_subelements[4]
    assert (
        qpointset_xml.tag == "qpointset"
    ), "Fifth xs subelement should be qpointset"
    assert qpointset_xml.items() == []
    qpoints = list(qpointset_xml)
    assert len(qpoints) == 1, "Excpect 1 qpoint"
    assert qpoints[0].tag == "qpoint"
    assert qpoints[0].items() == []
    assert qpoints[0].text == "0.0 0.0 0.0"


def test_exciting_own_parameters(tmp_path: Path, mp_Structure_mp390: Structure):
    own_cards = {
        "structure": {},
        "groundstate": {
            "nempty": "42",
        },
        "xs": {
            "xstype": "BSE",
            "nempty": "43",
            "energywindow": {"points": "44"},
            "BSE": {},
        },
    }

    exciting_parameters = EXCITINGParameters(
        own_cards, "keks", ["piep", "blub"], name="Test"
    )
    assert exciting_parameters.name == "Test"

    write_results = exciting_parameters.write(
        tmp_path, structure_uc=mp_Structure_mp390, sites=[1, 2, 3]
    )
    assert write_results == {"pass": True, "errors": {}}

    o1_path = tmp_path / "002_O"
    assert o1_path.exists(), "Directory for absorbing site 2 was not generated."
    input_file_o1 = o1_path / "input.xml"
    assert input_file_o1.exists(), "Input file not found."

    o2_path = tmp_path / "003_O"
    assert o2_path.exists(), "Directory for absorbing site 3 was not generated."
    input_file_o2 = o2_path / "input.xml"
    assert input_file_o2.exists(), "Input file not found."

    ti_path = tmp_path / "001_Ti"
    assert ti_path.exists(), "Directory for absorbing site 1 was not generated."
    input_file_ti = ti_path / "input.xml"
    assert input_file_ti.exists(), "Input file not found."

    root = ElementTree.parse(input_file_ti).getroot()

    assert root.tag == "input"

    subelements = list(root)
    assert len(subelements) == 4

    structure_xml = subelements[1]
    assert structure_xml.tag == "structure"
    assert (
        structure_xml.get("speciespath") == "keks"
    ), "species path set to keks"

    groundstate_xml = subelements[2]
    assert groundstate_xml.tag == "groundstate"
    assert (
        len(groundstate_xml.keys()) == 2
    ), "Excpect groundstate to have 2 attributes."
    assert set(groundstate_xml.keys()) == {
        "ngridk",
        "nempty",
    }, "Excpected groundstate keys differ from found keys."
    assert groundstate_xml.get("ngridk") == "3 3 4"
    assert groundstate_xml.get("nempty") == "42"

    xs_xml = subelements[3]
    assert xs_xml.tag == "xs"
    assert len(xs_xml.keys()) == 5, "Excpect xs to have 4 attributes."
    assert set(xs_xml.keys()) == {
        "ngridk",
        "nempty",
        "ngridq",
        "xstype",
        "gqmax",
    }, "Excpected xs keys differ from found keys."
    assert xs_xml.get("gqmax") == "4.0"
    assert xs_xml.get("ngridk") == "3 3 4"
    assert xs_xml.get("xstype") == "BSE"
    assert xs_xml.get("nempty") == "43"
    assert xs_xml.get("ngridq") == "3 3 4"

    xs_subelements = list(xs_xml)
    assert len(xs_subelements) == 3, "Expect xs tree to have 3 sub-elements"

    energywindow_xml = xs_subelements[0]
    assert (
        energywindow_xml.tag == "energywindow"
    ), "First xs subelement should be energywindow"
    assert set(energywindow_xml.keys()) == {"points"}
    assert energywindow_xml.get("points") == "44"

    plan_xml = xs_subelements[1]
    assert plan_xml.tag == "plan", "Second xs subelement should be plan"
    assert plan_xml.items() == []
    doonlys = list(plan_xml)
    assert len(doonlys) == 2, "Excpect 2 plan elements"
    assert doonlys[0].tag == "doonly"
    assert doonlys[0].items() == [("task", "piep")]
    assert doonlys[1].tag == "doonly"
    assert doonlys[1].items() == [("task", "blub")]

    bse_xml = xs_subelements[2]
    assert bse_xml.tag == "BSE", "Third xs subelement should be BSE"
    assert set(bse_xml.keys()) == {
        "nstlxas",
        "xasedge",
        "xasspecies",
        "xasatom",
    }
    assert bse_xml.get("nstlxas") == "1 59"
    assert bse_xml.get("xasspecies") == "1"
    assert bse_xml.get("xasatom") == "2"

    bse_xml_o1 = list(list(ElementTree.parse(input_file_o1).getroot())[3])[2]
    assert bse_xml_o1.get("xasspecies") == "1"
    assert bse_xml_o1.get("xasatom") == "3"

    bse_xml_o2 = list(list(ElementTree.parse(input_file_o2).getroot())[3])[2]
    assert bse_xml_o2.get("xasspecies") == "1"
    assert bse_xml_o2.get("xasatom") == "4"

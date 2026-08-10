"""Basic CLI regression tests for dmc.cli.epical."""

import pytest

from cpgmodule.cli.epical import build_parser, main


def test_build_parser():
    parser = build_parser()
    assert parser.prog


def test_top_level_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["-h"])
    assert excinfo.value.code == 0
    assert "Horvath13" in capsys.readouterr().out


def test_version(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert "epical" in capsys.readouterr().out


def test_horvath13_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["Horvath13", "-h"])
    assert excinfo.value.code == 0
    output = capsys.readouterr().out
    assert "Input_file" in output
    assert "--impute" in output
    assert "--log" in output


def test_dunedinpace_help_has_no_log(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["DunedinPACE", "-h"])
    assert excinfo.value.code == 0
    output = capsys.readouterr().out
    assert "Input_file" in output
    assert "--impute" in output
    assert "--log" not in output


def test_mouse_clock_has_genome_option():
    args = build_parser().parse_args(["WLMT", "input.tsv", "--genome", "mm39"])
    assert args.command == "WLMT"
    assert args.genome == "mm39"


def test_mammalian_clock_has_species_option():
    args = build_parser().parse_args(
        ["mammClock2", "input.tsv", "--species", "mouse"]
    )
    assert args.command == "mammClock2"
    assert args.species == "mouse"


def test_epm_specific_arguments():
    args = build_parser().parse_args(
        ["EPM", "beta.tsv", "meta.tsv", "--pcc", "0.9", "--kfold", "5"]
    )
    assert args.command == "EPM"
    assert args.input == "beta.tsv"
    assert args.meta == "meta.tsv"
    assert args.pcc == pytest.approx(0.9)
    assert args.kfold == 5

"""Tests for cpgmodule.cli.beta_impute."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cpgmodule.cli import beta_impute


def make_input(path: Path) -> Path:
    df = pd.DataFrame(
        {
            "S1": [0.10, np.nan, 0.30],
            "S2": [0.20, 0.40, np.nan],
            "S3": [0.30, 0.50, 0.70],
        },
        index=["cg1", "cg2", "cg3"],
    )
    df.to_csv(path, sep="\t", na_rep="NA")
    return path


def read_output(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def get_subcommands(parser: argparse.ArgumentParser) -> set[str]:
    action = next(
        a
        for a in parser._actions
        if isinstance(a, argparse._SubParsersAction)
    )
    return set(action.choices)


def test_build_parser():
    parser = beta_impute.build_parser()
    assert isinstance(parser, argparse.ArgumentParser)


def test_expected_subcommands():
    parser = beta_impute.build_parser()

    expected = {
        "toy",
        "insertna",
        "dropna",
        "countna",
        "constant",
        "mean",
        "median",
        "min",
        "max",
        "rand",
        "refknn",
        "mw",
        "knn",
        "buck",
        "rf",
        "softimpute",
        "morel",
        "gnn",
    }

    assert get_subcommands(parser) == expected


def test_top_level_help(capsys):
    parser = beta_impute.build_parser()

    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["-h"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "Impute missing values" in captured.out
    assert "softimpute" in captured.out
    assert "morel" in captured.out
    assert "gnn" in captured.out


def test_version(capsys):
    parser = beta_impute.build_parser()

    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["--version"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "nafiller" in captured.out


@pytest.mark.parametrize(
    "command",
    [
        "toy",
        "insertna",
        "dropna",
        "countna",
        "constant",
        "mean",
        "median",
        "min",
        "max",
        "rand",
        "refknn",
        "mw",
        "knn",
        "buck",
        "rf",
        "softimpute",
        "morel",
        "gnn",
    ],
)
def test_subcommand_help(command, capsys):
    parser = beta_impute.build_parser()

    with pytest.raises(SystemExit) as exc:
        parser.parse_args([command, "-h"])

    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "usage:" in captured.out.lower()


def test_toy_command(tmp_path):
    outfile = tmp_path / "toy.tsv"

    rc = beta_impute.main([
        "toy",
        str(outfile),
        "--rows", "6",
        "--cols", "4",
        "--missingness", "0.25",
        "--seed", "123",
    ])

    assert rc == 0
    assert outfile.exists()

    df = read_output(outfile)
    assert df.shape == (6, 4)
    assert int(df.isna().sum().sum()) > 0


def test_constant_imputation(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    outfile = tmp_path / "constant.tsv"

    rc = beta_impute.main([
        "constant",
        str(infile),
        str(outfile),
        "0.25",
    ])

    assert rc == 0
    df = read_output(outfile)

    assert df.loc["cg2", "S1"] == pytest.approx(0.25)
    assert df.loc["cg3", "S2"] == pytest.approx(0.25)
    assert int(df.isna().sum().sum()) == 0


def test_mean_imputation_by_row(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    outfile = tmp_path / "mean.tsv"

    rc = beta_impute.main([
        "mean",
        str(infile),
        str(outfile),
        "--axis", "index",
    ])

    assert rc == 0
    df = read_output(outfile)

    # cg2 mean of observed values: (0.40 + 0.50) / 2 = 0.45
    assert df.loc["cg2", "S1"] == pytest.approx(0.45)

    # cg3 mean of observed values: (0.30 + 0.70) / 2 = 0.50
    assert df.loc["cg3", "S2"] == pytest.approx(0.50)


def test_dropna_rows(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    outfile = tmp_path / "dropna.tsv"

    rc = beta_impute.main([
        "dropna",
        str(infile),
        str(outfile),
        "--axis", "rows",
    ])

    assert rc == 0
    df = read_output(outfile)

    assert list(df.index) == ["cg1"]


def test_countna_reports(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    row_report = tmp_path / "rows.tsv"
    col_report = tmp_path / "cols.tsv"

    rc = beta_impute.main([
        "countna",
        str(infile),
        "--row-report", str(row_report),
        "--column-report", str(col_report),
    ])

    assert rc == 0
    assert row_report.exists()
    assert col_report.exists()


def test_overwrite_protection(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    outfile = tmp_path / "constant.tsv"

    first = beta_impute.main([
        "constant",
        str(infile),
        str(outfile),
        "0.1",
    ])
    second = beta_impute.main([
        "constant",
        str(infile),
        str(outfile),
        "0.2",
    ])

    assert first == 0
    assert second == 1


def test_overwrite_allows_replacement(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    outfile = tmp_path / "constant.tsv"

    assert beta_impute.main([
        "constant",
        str(infile),
        str(outfile),
        "0.1",
    ]) == 0

    assert beta_impute.main([
        "constant",
        str(infile),
        str(outfile),
        "0.2",
        "--overwrite",
    ]) == 0

    df = read_output(outfile)
    assert df.loc["cg2", "S1"] == pytest.approx(0.2)


def test_truth_evaluation_path(tmp_path):
    infile = make_input(tmp_path / "input.tsv")
    outfile = tmp_path / "constant.tsv"
    truth = tmp_path / "truth.tsv"

    truth_df = pd.DataFrame(
        {
            "S1": [0.10, 0.45, 0.30],
            "S2": [0.20, 0.40, 0.55],
            "S3": [0.30, 0.50, 0.70],
        },
        index=["cg1", "cg2", "cg3"],
    )
    truth_df.to_csv(truth, sep="\t")

    rc = beta_impute.main([
        "constant",
        str(infile),
        str(outfile),
        "0.5",
        "--truth", str(truth),
    ])

    assert rc == 0
    assert outfile.exists()


def test_load_morel_groups(tmp_path):
    group_file = tmp_path / "groups.json"
    group_file.write_text(
        '{"group_1": ["S1", "S2"], "group_2": ["S3", "S4"]}',
        encoding="utf-8",
    )

    groups = beta_impute.load_morel_groups(group_file)

    assert groups == {
        "group_1": ["S1", "S2"],
        "group_2": ["S3", "S4"],
    }


def test_load_morel_groups_invalid_json(tmp_path):
    group_file = tmp_path / "groups.json"
    group_file.write_text("{bad json", encoding="utf-8")

    with pytest.raises(ValueError, match="Invalid JSON"):
        beta_impute.load_morel_groups(group_file)


def test_read_df_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        beta_impute.read_df(tmp_path / "missing.tsv")


def test_write_df_rejects_negative_decimal(tmp_path):
    df = pd.DataFrame({"S1": [0.1]}, index=["cg1"])

    with pytest.raises(ValueError, match="non-negative"):
        beta_impute.write_df(
            df,
            tmp_path / "out.tsv",
            decimal=-1,
        )

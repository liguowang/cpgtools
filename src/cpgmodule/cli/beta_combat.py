#!/usr/bin/env python3

# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: wangliguo78@gmail.com
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

"""Correct batch effects in DNA methylation beta-value matrices with ComBat.

The input beta matrix must be tab-delimited, with CpG identifiers in the first
column and sample identifiers in the header. The batch file must contain two
columns: sample identifier and batch/group identifier.

Missing beta values are temporarily imputed with K-nearest-neighbor imputation
before ComBat is applied. When missing values are present, the program writes
both a fully imputed, batch-corrected matrix and a second matrix in which the
original missing-value positions are restored.

Description
-----------
This program corrects batch effect using the "combat" algorithm:

W. Evan Johnson, et al, Adjusting batch effects in microarray expression data using empirical Bayes methods, Biostatistics, 2007.

Example of input data file
---------------------------
CpG_ID    Sample_01    Sample_02    Sample_03    Sample_04
cg_001    0.831035    0.878022    0.794427    0.880911
cg_002    0.249544    0.209949    0.234294    0.236680
cg_003    0.845065    0.843957    0.840184    0.824286
...

Example of batch file
-------------------------------
Sample_01,plate_1
Sample_02,plate_1
Sample_03,plate_2
Sample_04,plate_2
...

"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from combat.pycombat import pycombat
from sklearn.impute import KNNImputer

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class CombatError(RuntimeError):
    """Raised when a ComBat input or correction step is invalid."""


def read_beta_matrix(path: Path) -> pd.DataFrame:
    """Read and validate a tab-delimited beta-value matrix."""
    if not path.is_file():
        raise CombatError(f"Input beta matrix does not exist: {path}")

    try:
        data = pd.read_csv(path, sep="\t", index_col=0, compression="infer")
    except Exception as exc:
        raise CombatError(f"Cannot read beta matrix {path}: {exc}") from exc

    if data.empty:
        raise CombatError(f"Input beta matrix is empty: {path}")
    if data.shape[1] < 2:
        raise CombatError("ComBat requires at least two samples.")
    if not data.index.is_unique:
        raise CombatError("CpG identifiers in the beta matrix must be unique.")
    if not data.columns.is_unique:
        raise CombatError("Sample identifiers in the beta matrix must be unique.")

    data.columns = data.columns.map(str)
    data.index = data.index.map(str)
    data = data.apply(pd.to_numeric, errors="coerce")

    if data.notna().sum().sum() == 0:
        raise CombatError("The beta matrix contains no numeric values.")

    finite_values = data.to_numpy(dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size and ((finite_values < 0).any() or (finite_values > 1).any()):
        printlog(
            "Warning: beta values outside [0, 1] were detected; "
            "they will be retained during correction."
        )

    return data


def read_batch_file(path: Path) -> pd.Series:
    """Read a two-column sample-to-batch mapping file.

    The file may be comma-, tab-, or whitespace-delimited and may optionally
    contain a header row. Reading with ``header=None`` prevents the first sample
    assignment from being mistaken for a header in headerless files.
    """
    if not path.is_file():
        raise CombatError(f"Batch file does not exist: {path}")

    try:
        batch_table = pd.read_csv(
            path,
            sep=None,
            engine="python",
            comment="#",
            header=None,
            dtype=str,
        )
    except Exception as exc:
        raise CombatError(f"Cannot read batch file {path}: {exc}") from exc

    if batch_table.shape[1] < 2:
        raise CombatError("The batch file must contain at least two columns.")

    batch_table = batch_table.iloc[:, :2].copy()
    batch_table.columns = ["sample", "batch"]

    # Normalize identifiers and remove blank records. lstrip removes a possible
    # UTF-8 byte-order mark from the first field.
    batch_table["sample"] = (
        batch_table["sample"].fillna("").astype(str).str.strip().str.lstrip("\ufeff")
    )
    batch_table["batch"] = batch_table["batch"].fillna("").astype(str).str.strip()
    batch_table = batch_table.loc[
        (batch_table["sample"] != "") & (batch_table["batch"] != "")
    ].copy()

    # Remove an optional header row.
    if not batch_table.empty:
        first_sample = batch_table.iloc[0, 0].lower()
        first_batch = batch_table.iloc[0, 1].lower()
        if first_sample in {"sample", "sample_id", "sampleid"} and first_batch in {
            "group",
            "group_id",
            "groupid",
            "batch",
            "batch_id",
            "batchid",
        }:
            batch_table = batch_table.iloc[1:].copy()

    if batch_table.empty:
        raise CombatError("The batch file contains no sample assignments.")

    if batch_table["sample"].duplicated().any():
        duplicates = batch_table.loc[
            batch_table["sample"].duplicated(keep=False), "sample"
        ].unique()
        raise CombatError(
            "Sample identifiers in the batch file must be unique: "
            + ", ".join(map(str, duplicates[:10]))
        )

    return batch_table.set_index("sample")["batch"]


def align_batches(samples: Sequence[str], batches: pd.Series) -> pd.Series:
    """Align batch labels to beta-matrix sample order and validate membership."""
    sample_index = pd.Index(map(str, samples), name="sample")
    missing = sample_index.difference(batches.index)
    extra = batches.index.difference(sample_index)

    if len(missing):
        raise CombatError(
            "Samples missing from the batch file: " + ", ".join(missing[:10])
        )
    if len(extra):
        printlog(
            "Warning: batch assignments for samples absent from the beta matrix "
            "will be ignored: " + ", ".join(extra[:10])
        )

    aligned = batches.reindex(sample_index)
    if aligned.nunique() < 2:
        raise CombatError("ComBat requires at least two distinct batch groups.")

    batch_counts = aligned.value_counts()
    singleton_batches = batch_counts[batch_counts < 2]
    if not singleton_batches.empty:
        printlog(
            "Warning: batch groups with fewer than two samples were detected: "
            + ", ".join(
                f"{batch}={count}" for batch, count in singleton_batches.items()
            )
        )

    return aligned


def impute_missing_values(
    data: pd.DataFrame,
    n_neighbors: int,
    axis: int,
) -> pd.DataFrame:
    """Impute missing beta values with K-nearest-neighbor imputation."""
    if n_neighbors < 1:
        raise CombatError("--neighbors must be at least 1.")
    if axis not in (0, 1):
        raise CombatError("--axis must be 0 or 1.")

    working = data.T if axis == 1 else data
    max_neighbors = working.shape[0]
    if n_neighbors > max_neighbors:
        raise CombatError(
            f"--neighbors ({n_neighbors}) exceeds the number of available "
            f"observations ({max_neighbors}) for axis={axis}."
        )

    imputer = KNNImputer(n_neighbors=n_neighbors)
    imputed = imputer.fit_transform(working)

    # KNNImputer can drop all-missing features in older scikit-learn versions.
    if imputed.shape != working.shape:
        raise CombatError(
            "KNN imputation changed the matrix dimensions. This usually means "
            "that one or more rows or columns contain only missing values."
        )

    result = pd.DataFrame(imputed, index=working.index, columns=working.columns)
    return result.T if axis == 1 else result


def pick_colors(batch_labels: pd.Series) -> list[str]:
    """Return one plotting color per sample, preserving sample order."""
    unique_batches = list(dict.fromkeys(batch_labels.tolist()))
    available = list(mcolors.TABLEAU_COLORS.values()) + list(mcolors.CSS4_COLORS.values())
    if len(unique_batches) > len(available):
        raise CombatError(
            f"Too many batch groups for plotting: {len(unique_batches)}."
        )

    color_by_batch = dict(zip(unique_batches, available))
    return [color_by_batch[label] for label in batch_labels]


def box_plot(
    data: pd.DataFrame,
    sample_colors: Sequence[str],
    output_path: Path,
    ylabel: str = "Beta values",
    title: str = "",
) -> None:
    """Write a sample-wise beta-value box plot."""
    figure_width = max(10.0, min(30.0, 0.35 * data.shape[1] + 4.0))
    figure, axis = plt.subplots(figsize=(figure_width, 6.0))
    boxplot = axis.boxplot(
        [data[column].dropna().to_numpy() for column in data.columns],
        patch_artist=True,
        tick_labels=data.columns,
        showfliers=False,
    )

    for patch, color in zip(boxplot["boxes"], sample_colors):
        patch.set_facecolor(color)

    axis.set_xticklabels(data.columns, rotation=90, fontsize=8)
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300)
    plt.close(figure)


def run_combat(
    beta_data: pd.DataFrame,
    batch_labels: pd.Series,
) -> pd.DataFrame:
    """Apply ComBat and return a DataFrame with preserved labels."""
    try:
        corrected = pycombat(beta_data, batch_labels.tolist())
    except Exception as exc:
        raise CombatError(f"ComBat correction failed: {exc}") from exc

    if isinstance(corrected, pd.DataFrame):
        result = corrected.copy()
        result.index = beta_data.index
        result.columns = beta_data.columns
    else:
        corrected_array = np.asarray(corrected, dtype=float)
        if corrected_array.shape != beta_data.shape:
            raise CombatError(
                "ComBat returned a matrix with unexpected dimensions: "
                f"{corrected_array.shape}; expected {beta_data.shape}."
            )
        result = pd.DataFrame(
            corrected_array,
            index=beta_data.index,
            columns=beta_data.columns,
        )

    return result


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Correct batch effects in a DNA methylation beta-value matrix "
            "using ComBat."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Tab-delimited beta-value matrix with CpG IDs in column 1 and "
            "sample IDs in the header."
        ),
    )
    parser.add_argument(
        "-g",
        "--group",
        required=True,
        type=Path,
        dest="group_file",
        help=(
            "Two-column batch file containing sample ID and batch/group ID. "
            "Comma- or tab-delimited files are accepted."
        ),
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        "--output",
        required=True,
        type=Path,
        dest="out_prefix",
        help="Output filename prefix, optionally including a directory.",
    )
    parser.add_argument(
        "-k",
        "--neighbors",
        type=int,
        default=3,
        dest="n_neighbors",
        help="Number of neighbors used by KNN imputation.",
    )
    parser.add_argument(
        "--axis",
        type=int,
        choices=(0, 1),
        default=1,
        dest="axis_choice",
        help=(
            "KNN-imputation orientation: 1 finds neighboring samples after "
            "transposing the matrix; 0 finds neighboring CpGs."
        ),
    )
    parser.add_argument(
        "--no_plot",
        action="store_true",
        help="Do not generate box plots before and after ComBat correction.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    corrected_path = Path(f"{args.out_prefix}.combat.tsv")
    restored_na_path = Path(f"{args.out_prefix}.combat_withNAs.tsv")
    before_plot_path = Path(f"{args.out_prefix}.boxplot.png")
    after_plot_path = Path(f"{args.out_prefix}.boxplot_combat.png")

    try:
        printlog(f'Reading input file: "{args.input_file}" ...')
        original = read_beta_matrix(args.input_file)

        printlog(f'Reading batch file: "{args.group_file}" ...')
        batches = align_batches(original.columns, read_batch_file(args.group_file))
        batch_counts = batches.value_counts(sort=False)
        printlog(
            "Batch composition: "
            + ", ".join(f"{batch}={count}" for batch, count in batch_counts.items())
        )

        missing_mask = original.isna()
        total_missing = int(missing_mask.sum().sum())
        printlog(f"The input matrix contains {total_missing} missing values.")

        if total_missing:
            printlog(
                f"Imputing missing values with KNN (k={args.n_neighbors}, "
                f"axis={args.axis_choice}) ..."
            )
            working = impute_missing_values(
                original,
                n_neighbors=args.n_neighbors,
                axis=args.axis_choice,
            )
        else:
            working = original.copy()

        sample_colors = pick_colors(batches)
        if not args.no_plot:
            printlog(f'Writing pre-correction box plot: "{before_plot_path}" ...')
            box_plot(
                working,
                sample_colors=sample_colors,
                output_path=before_plot_path,
                title="Before batch-effect correction",
            )

        printlog("Removing batch effects with ComBat ...")
        corrected = run_combat(working, batches)

        if not args.no_plot:
            printlog(f'Writing post-correction box plot: "{after_plot_path}" ...')
            box_plot(
                corrected,
                sample_colors=sample_colors,
                output_path=after_plot_path,
                title="After batch-effect correction",
            )

        printlog(
            f'Writing corrected matrix with imputed values: "{corrected_path}" ...'
        )
        corrected.to_csv(corrected_path, sep="\t", float_format="%.6f")

        if total_missing:
            corrected_with_na = corrected.mask(missing_mask)
            printlog(
                "Writing corrected matrix with original missing values restored: "
                f'"{restored_na_path}" ...'
            )
            corrected_with_na.to_csv(
                restored_na_path,
                sep="\t",
                na_rep="NA",
                float_format="%.6f",
            )

    except CombatError as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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

"""Select the top N CpGs according to a row-wise summary score.

The input matrix must contain CpG identifiers in the first column and sample
values in the remaining columns. CpGs can be ranked by standard deviation,
variance, mean, median, or median absolute deviation (MAD).

Two output files are written:

* ``PREFIX.ranked.tsv``: all retained CpGs with their ranking score.
* ``PREFIX.topN.tsv``: the top-ranked CpGs without the score column.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class TopNError(RuntimeError):
    """Raised when input validation or ranking fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_topN",
        description=(
            "Rank CpGs by a row-wise score and retain the top N features for "
            "downstream analyses such as PCA or clustering."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input tabular matrix. The first column contains CpG IDs and "
            "remaining columns contain sample values. Compressed input is supported."
        ),
    )
    parser.add_argument(
        "-c",
        "--count",
        "--top",
        dest="top_n",
        type=int,
        default=1000,
        help="Number of top-ranked CpGs to retain.",
    )
    parser.add_argument(
        "-s",
        "--score",
        choices=("std", "var", "mean", "median", "mad"),
        default="std",
        help="Row-wise score used to rank CpGs.",
    )
    parser.add_argument(
        "--ascending",
        action="store_true",
        help="Sort from smallest to largest score instead of descending order.",
    )
    parser.add_argument(
        "--na_policy",
        choices=("drop", "omit"),
        default="drop",
        help=(
            "How to handle missing values. 'drop' removes any CpG containing a "
            "missing value; 'omit' computes the score from available values."
        ),
    )
    parser.add_argument(
        "--min_valid",
        type=int,
        default=2,
        help=(
            "Minimum number of non-missing sample values required per CpG when "
            "--na_policy omit is used."
        ),
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        "--output",
        dest="out_prefix",
        required=True,
        type=Path,
        help="Output filename prefix, optionally including a directory.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Validate command-line arguments."""
    if not args.input_file.is_file():
        parser.error(f"Input file does not exist: {args.input_file}")
    if args.top_n < 1:
        parser.error("--count/--top must be at least 1")
    if args.min_valid < 1:
        parser.error("--min_valid must be at least 1")


def read_matrix(path: Path) -> pd.DataFrame:
    """Read and validate the input feature matrix."""
    printlog(f'Reading input matrix: "{path}"')
    try:
        data = pd.read_csv(
            path,
            sep=None,
            engine="python",
            compression="infer",
            index_col=0,
        )
    except Exception as exc:
        raise TopNError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] < 1:
        raise TopNError("Input must contain at least one sample column.")

    if data.index.has_duplicates:
        duplicates = data.index[data.index.duplicated()].astype(str).unique().tolist()
        raise TopNError(
            "Duplicate CpG identifiers were found: " + ", ".join(duplicates[:10])
        )

    if data.columns.duplicated().any():
        duplicates = data.columns[data.columns.duplicated()].astype(str).tolist()
        raise TopNError(
            "Duplicate sample identifiers were found: " + ", ".join(duplicates[:10])
        )

    numeric = data.apply(pd.to_numeric, errors="coerce")
    if numeric.notna().sum().sum() == 0:
        raise TopNError("No numeric sample values were found in the input matrix.")

    return numeric


def median_absolute_deviation(data: pd.DataFrame, skipna: bool) -> pd.Series:
    """Calculate row-wise median absolute deviation."""
    medians = data.median(axis=1, skipna=skipna)
    deviations = data.sub(medians, axis=0).abs()
    return deviations.median(axis=1, skipna=skipna)


def compute_scores(
    data: pd.DataFrame,
    score_type: str,
    na_policy: str,
    min_valid: int,
) -> tuple[pd.DataFrame, pd.Series]:
    """Filter CpGs and calculate row-wise ranking scores."""
    original_count = len(data)

    if na_policy == "drop":
        filtered = data.dropna(axis=0, how="any")
        removed = original_count - len(filtered)
        printlog(f"Removed {removed} CpGs containing missing values.")
        skipna = False
    else:
        valid_counts = data.notna().sum(axis=1)
        filtered = data.loc[valid_counts >= min_valid]
        removed = original_count - len(filtered)
        printlog(
            f"Removed {removed} CpGs with fewer than {min_valid} valid values."
        )
        skipna = True

    if filtered.empty:
        raise TopNError("No CpGs remain after missing-value filtering.")

    if score_type == "std":
        scores = filtered.std(axis=1, skipna=skipna, ddof=1)
    elif score_type == "var":
        scores = filtered.var(axis=1, skipna=skipna, ddof=1)
    elif score_type == "mean":
        scores = filtered.mean(axis=1, skipna=skipna)
    elif score_type == "median":
        scores = filtered.median(axis=1, skipna=skipna)
    elif score_type == "mad":
        scores = median_absolute_deviation(filtered, skipna=skipna)
    else:  # Defensive check; argparse already restricts this value.
        raise TopNError(f"Unsupported score type: {score_type}")

    finite_scores = np.isfinite(scores.to_numpy(dtype=float))
    if not finite_scores.all():
        bad_count = int((~finite_scores).sum())
        printlog(f"Removed {bad_count} CpGs with undefined ranking scores.")
        filtered = filtered.loc[finite_scores]
        scores = scores.loc[finite_scores]

    if filtered.empty:
        raise TopNError("No CpGs have finite ranking scores.")

    return filtered, scores.astype(float)


def rank_features(
    data: pd.DataFrame,
    scores: pd.Series,
    ascending: bool,
) -> pd.DataFrame:
    """Return the feature matrix sorted by score with score and rank columns."""
    ranked = data.copy()
    ranked.insert(0, "Score", scores)
    ranked = ranked.sort_values(
        by="Score",
        ascending=ascending,
        kind="mergesort",
    )
    ranked.insert(1, "Rank", np.arange(1, len(ranked) + 1, dtype=int))
    return ranked


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        matrix = read_matrix(args.input_file)
        filtered, scores = compute_scores(
            data=matrix,
            score_type=args.score,
            na_policy=args.na_policy,
            min_valid=args.min_valid,
        )
        ranked = rank_features(filtered, scores, ascending=args.ascending)

        selected_count = min(args.top_n, len(ranked))
        if args.top_n > len(ranked):
            printlog(
                f"Requested {args.top_n} CpGs, but only {len(ranked)} are available; "
                f"retaining all {len(ranked)} CpGs."
            )

        ranked_path = Path(f"{args.out_prefix}.ranked.tsv")
        top_path = Path(f"{args.out_prefix}.topN.tsv")
        ids_path = Path(f"{args.out_prefix}.topN_ids.txt")

        printlog(f'Writing ranked feature table: "{ranked_path}"')
        ranked.to_csv(ranked_path, sep="\t", float_format="%.6f")

        top_matrix = ranked.iloc[:selected_count].drop(columns=["Score", "Rank"])
        printlog(f'Writing top {selected_count} feature matrix: "{top_path}"')
        top_matrix.to_csv(top_path, sep="\t", float_format="%.6f")

        printlog(f'Writing selected CpG identifiers: "{ids_path}"')
        ids_path.write_text(
            "\n".join(map(str, top_matrix.index)) + "\n",
            encoding="utf-8",
        )

    except TopNError as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

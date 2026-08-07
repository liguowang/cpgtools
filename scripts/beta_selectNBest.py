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

"""Select the top-ranked CpG features for a classification phenotype.

Supported scoring methods
-------------------------
* ANOVA F statistic (``anova``)
* Mutual information (``mi``)
* Chi-square statistic (``chisq``; requires non-negative values)

The input methylation matrix must contain CpG identifiers in the first column
and sample identifiers in the remaining columns. The group file must contain
sample identifiers in the first column and class labels in the second column.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, chi2, f_classif, mutual_info_classif

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class FeatureSelectionError(RuntimeError):
    """Raised when feature-selection inputs or results are invalid."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_selectNBest",
        description="Select top CpG features for a categorical phenotype.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Tab-delimited beta-value matrix. The first column contains CpG IDs "
            "and the remaining columns contain samples. Compressed input is supported."
        ),
    )
    parser.add_argument(
        "-g",
        "--group",
        required=True,
        type=Path,
        help=(
            "Two-column sample-group file. The first column contains sample IDs "
            "and the second contains categorical labels. CSV and TSV are supported."
        ),
    )
    parser.add_argument(
        "-k",
        "-c",
        "--top",
        "--topK",
        dest="top_k",
        type=int,
        default=100,
        help="Number of top-ranked CpG features to select.",
    )
    parser.add_argument(
        "-s",
        "--score_function",
        "--score-function",
        dest="score_function",
        choices=("anova", "mi", "chisq"),
        default="chisq",
        help="Feature-scoring method.",
    )
    parser.add_argument(
        "--random_state",
        type=int,
        default=0,
        help="Random seed used by mutual-information scoring.",
    )
    parser.add_argument(
        "--keep_constant",
        action="store_true",
        help="Do not remove zero-variance CpGs before feature selection.",
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


def read_beta_matrix(path: Path) -> pd.DataFrame:
    """Read and validate a methylation matrix."""
    if not path.is_file():
        raise FeatureSelectionError(f"Input methylation file does not exist: {path}")

    try:
        data = pd.read_csv(path, sep="\t", compression="infer")
    except Exception as exc:
        raise FeatureSelectionError(f"Cannot read methylation matrix {path}: {exc}") from exc

    if data.shape[1] < 3:
        raise FeatureSelectionError(
            "The methylation matrix must contain one CpG column and at least two samples."
        )

    cpg_column = data.columns[0]
    data[cpg_column] = data[cpg_column].astype(str)

    if data[cpg_column].duplicated().any():
        duplicates = data.loc[data[cpg_column].duplicated(), cpg_column].head(5).tolist()
        raise FeatureSelectionError(
            "Duplicate CpG identifiers were found, for example: " + ", ".join(duplicates)
        )

    sample_names = [str(column) for column in data.columns[1:]]
    if len(sample_names) != len(set(sample_names)):
        raise FeatureSelectionError("Sample identifiers in the methylation matrix are not unique.")

    values = data.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    values.index = data[cpg_column]
    values.index.name = "CpG_ID"
    values.columns = sample_names
    return values


def read_group_file(path: Path) -> pd.Series:
    """Read a two-column sample-group file with or without a header."""
    if not path.is_file():
        raise FeatureSelectionError(f"Group file does not exist: {path}")

    try:
        group = pd.read_csv(path, sep=None, engine="python", header=None, dtype=str)
    except Exception as exc:
        raise FeatureSelectionError(f"Cannot read group file {path}: {exc}") from exc

    if group.shape[1] < 2:
        raise FeatureSelectionError("The group file must contain at least two columns.")

    group = group.iloc[:, :2].copy()
    group.columns = ["Sample_ID", "Group_ID"]
    group["Sample_ID"] = group["Sample_ID"].str.strip().str.lstrip("\ufeff")
    group["Group_ID"] = group["Group_ID"].str.strip()

    if not group.empty:
        first_sample = group.iloc[0, 0].lower()
        first_group = group.iloc[0, 1].lower()
        if first_sample in {"sample", "sample_id", "sampleid"} and first_group in {
            "group",
            "group_id",
            "groupid",
            "class",
            "label",
        }:
            group = group.iloc[1:].copy()

    group = group.dropna(subset=["Sample_ID", "Group_ID"])
    group = group[(group["Sample_ID"] != "") & (group["Group_ID"] != "")]

    if group.empty:
        raise FeatureSelectionError("No valid sample-group assignments were found.")

    if group["Sample_ID"].duplicated().any():
        duplicates = group.loc[group["Sample_ID"].duplicated(), "Sample_ID"].tolist()
        raise FeatureSelectionError(
            "Duplicate sample identifiers in group file: " + ", ".join(duplicates[:10])
        )

    labels = group.set_index("Sample_ID")["Group_ID"]
    labels.index = labels.index.astype(str)
    labels.name = "Group_ID"
    return labels


def prepare_data(
    beta: pd.DataFrame,
    labels: pd.Series,
    keep_constant: bool,
) -> tuple[pd.DataFrame, pd.Series, int, int]:
    """Align samples, remove incomplete CpGs, and optionally drop constants."""
    missing = [sample for sample in beta.columns if sample not in labels.index]
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise FeatureSelectionError(f"Samples missing from the group file: {preview}{suffix}")

    labels = labels.loc[beta.columns]
    if labels.nunique() < 2:
        raise FeatureSelectionError("At least two distinct groups are required.")

    complete = beta.dropna(axis=0, how="any")
    removed_missing = beta.shape[0] - complete.shape[0]

    if complete.empty:
        raise FeatureSelectionError("No complete CpG rows remain after removing missing values.")

    removed_constant = 0
    if not keep_constant:
        variances = complete.var(axis=1, ddof=0)
        keep = np.isfinite(variances) & (variances > 0)
        removed_constant = int((~keep).sum())
        complete = complete.loc[keep]

    if complete.empty:
        raise FeatureSelectionError("No variable CpG features remain for selection.")

    # Samples x features, as expected by scikit-learn.
    return complete.T, labels, removed_missing, removed_constant


def build_selector(method: str, top_k: int, random_state: int) -> SelectKBest:
    """Construct the requested scikit-learn feature selector."""
    if method == "anova":
        return SelectKBest(score_func=f_classif, k=top_k)
    if method == "mi":
        score_func = lambda x, y: mutual_info_classif(
            x,
            y,
            random_state=random_state,
        )
        return SelectKBest(score_func=score_func, k=top_k)
    return SelectKBest(score_func=chi2, k=top_k)


def run_feature_selection(
    matrix: pd.DataFrame,
    labels: pd.Series,
    method: str,
    top_k: int,
    random_state: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run feature selection and return selected data and ranked statistics."""
    if top_k < 1:
        raise FeatureSelectionError("--top must be at least 1.")

    number_features = matrix.shape[1]
    selected_k = min(top_k, number_features)

    if method == "chisq" and (matrix.to_numpy(dtype=float) < 0).any():
        raise FeatureSelectionError(
            "Chi-square scoring requires non-negative values. Use beta values or "
            "choose --score_function anova/mi for M values."
        )

    selector = build_selector(method, selected_k, random_state)
    selector.fit(matrix.to_numpy(dtype=float), labels.to_numpy())

    scores = np.asarray(selector.scores_, dtype=float)
    pvalues = (
        np.asarray(selector.pvalues_, dtype=float)
        if getattr(selector, "pvalues_", None) is not None
        else np.full(number_features, np.nan)
    )
    selected_mask = selector.get_support()

    stats = pd.DataFrame(
        {
            "CpG_ID": matrix.columns,
            "score": scores,
            "p_value": pvalues,
            "selected": selected_mask,
        }
    )
    stats["rank"] = stats["score"].rank(
        method="min",
        ascending=False,
        na_option="bottom",
    ).astype("Int64")
    stats = stats.sort_values(
        ["selected", "score", "CpG_ID"],
        ascending=[False, False, True],
        na_position="last",
    ).reset_index(drop=True)

    selected_features = stats.loc[stats["selected"], "CpG_ID"].tolist()
    selected_data = matrix.loc[:, selected_features].T
    selected_data.index.name = "CpG_ID"

    return selected_data, stats


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        printlog(f'Reading input file: "{args.input_file}"')
        beta = read_beta_matrix(args.input_file)

        printlog(f'Reading group file: "{args.group}"')
        labels = read_group_file(args.group)

        matrix, labels, removed_missing, removed_constant = prepare_data(
            beta,
            labels,
            keep_constant=args.keep_constant,
        )

        printlog(f"Removed {removed_missing} CpGs containing missing values.")
        if not args.keep_constant:
            printlog(f"Removed {removed_constant} zero-variance CpGs.")
        printlog(
            f"Using {matrix.shape[0]} samples and {matrix.shape[1]} CpG features "
            f"across {labels.nunique()} groups."
        )

        selected_data, stats = run_feature_selection(
            matrix=matrix,
            labels=labels,
            method=args.score_function,
            top_k=args.top_k,
            random_state=args.random_state,
        )

        args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
        selected_path = Path(f"{args.out_prefix}.selected_features.tsv")
        scores_path = Path(f"{args.out_prefix}.feature_scores.tsv")
        cpg_path = Path(f"{args.out_prefix}.selected_cpgs.txt")

        printlog(f'Writing selected feature matrix: "{selected_path}"')
        selected_data.to_csv(
            selected_path,
            sep="\t",
            index=True,
            float_format="%.6f",
        )

        printlog(f'Writing ranked feature scores: "{scores_path}"')
        stats.to_csv(
            scores_path,
            sep="\t",
            index=False,
            float_format="%.10g",
        )

        printlog(f'Writing selected CpG list: "{cpg_path}"')
        stats.loc[stats["selected"], "CpG_ID"].to_csv(
            cpg_path,
            index=False,
            header=False,
        )

        printlog(f"Selected {selected_data.shape[0]} CpG features.")

    except FeatureSelectionError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: unexpected failure: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

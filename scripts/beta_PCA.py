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

"""Perform principal component analysis on DNA methylation samples.

The input matrix must contain CpG identifiers in the first column and sample
beta values in subsequent columns. Samples are represented by rows after the
matrix is transposed. CpGs containing missing values across the selected
samples are removed before PCA.

A group file may be supplied to color samples in the PCA plot. It must contain
two columns: sample ID and group ID. Comma- and tab-delimited files, with or
without a header, are supported.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class PCAError(RuntimeError):
    """Raised when PCA input, metadata, or model fitting is invalid."""


def read_beta_matrix(path: Path) -> pd.DataFrame:
    """Read and validate a DNA methylation beta-value matrix."""
    if not path.is_file():
        raise PCAError(f"Input file does not exist: {path}")

    try:
        data = pd.read_csv(path, sep=None, engine="python", compression="infer")
    except Exception as exc:
        raise PCAError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] < 3:
        raise PCAError(
            "Input must contain a CpG column and at least two sample columns."
        )

    data = data.rename(columns={data.columns[0]: "CpG_ID"})
    data["CpG_ID"] = data["CpG_ID"].astype(str).str.strip()

    if data["CpG_ID"].duplicated().any():
        duplicate_count = int(data["CpG_ID"].duplicated().sum())
        printlog(
            f"Warning: removing {duplicate_count} duplicate CpG identifiers."
        )
        data = data.drop_duplicates(subset="CpG_ID", keep="first")

    sample_names = [str(column).strip() for column in data.columns[1:]]
    if len(sample_names) != len(set(sample_names)):
        raise PCAError("Sample IDs in the beta matrix are not unique.")
    data.columns = ["CpG_ID", *sample_names]

    for column in sample_names:
        data[column] = pd.to_numeric(data[column], errors="coerce")

    if data.iloc[:, 1:].notna().sum().sum() == 0:
        raise PCAError("No numeric beta values were found in the input matrix.")

    return data


def read_group_file(path: Path) -> pd.Series:
    """Read a two-column sample/group file with or without a header."""
    if not path.is_file():
        raise PCAError(f"Group file does not exist: {path}")

    try:
        group = pd.read_csv(
            path,
            sep=None,
            engine="python",
            header=None,
            dtype=str,
            comment="#",
        )
    except Exception as exc:
        raise PCAError(f"Cannot read group file {path}: {exc}") from exc

    group = group.dropna(how="all")
    if group.shape[1] < 2:
        raise PCAError("Group file must contain at least two columns.")

    group = group.iloc[:, :2].copy()
    group.columns = ["Sample_ID", "Group_ID"]
    for column in group.columns:
        group[column] = (
            group[column]
            .astype(str)
            .str.replace("\ufeff", "", regex=False)
            .str.strip()
        )

    if not group.empty:
        first_sample = group.iloc[0, 0].lower()
        first_group = group.iloc[0, 1].lower()
        sample_headers = {"sample", "sample_id", "sampleid", "id"}
        group_headers = {"group", "group_id", "groupid", "class", "batch"}
        if first_sample in sample_headers and first_group in group_headers:
            group = group.iloc[1:].copy()

    group = group[(group["Sample_ID"] != "") & (group["Group_ID"] != "")]
    if group.empty:
        raise PCAError("Group file contains no sample assignments.")

    if group["Sample_ID"].duplicated().any():
        duplicates = sorted(
            group.loc[group["Sample_ID"].duplicated(), "Sample_ID"].unique()
        )
        raise PCAError(
            "Sample IDs are not unique in the group file: "
            + ", ".join(duplicates[:10])
        )

    return group.set_index("Sample_ID")["Group_ID"]


def choose_colors(group_names: Sequence[str]) -> dict[str, object]:
    """Create a deterministic mapping from group labels to plot colors."""
    cmap = plt.get_cmap("tab20")
    count = len(group_names)
    if count == 1:
        return {group_names[0]: cmap(0)}
    return {
        name: cmap(index / max(1, count - 1))
        for index, name in enumerate(group_names)
    }


def plot_pca(
    scores: pd.DataFrame,
    explained_variance: np.ndarray,
    output_path: Path,
    label_samples: bool,
    marker: str,
    alpha: float,
    legend_location: str,
    width: float,
    height: float,
    dpi: int,
) -> None:
    """Write a two-dimensional PCA scatter plot."""
    figure, axis = plt.subplots(figsize=(width, height))

    group_names = list(dict.fromkeys(scores["Group_ID"].astype(str)))
    colors = choose_colors(group_names)

    for group_name in group_names:
        subset = scores[scores["Group_ID"] == group_name]
        axis.scatter(
            subset["PC1"],
            subset["PC2"],
            label=group_name,
            marker=marker,
            alpha=alpha,
            edgecolors="none" if marker != "o" else "black",
        )

        # Matplotlib assigns colors automatically per call; labels remain stable.
        if label_samples:
            for sample_id, row in subset.iterrows():
                axis.annotate(
                    str(sample_id),
                    (row["PC1"], row["PC2"]),
                    xytext=(3, -8),
                    textcoords="offset points",
                    fontsize=7,
                )

    axis.set_xlabel(f"PC1 ({explained_variance[0] * 100:.2f}% variance)")
    axis.set_ylabel(f"PC2 ({explained_variance[1] * 100:.2f}% variance)")
    axis.set_title("DNA methylation PCA")
    axis.legend(loc=legend_location, frameon=True)
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_PCA",
        description="Perform PCA on DNA methylation samples.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help="Input beta-value matrix; CpGs in column 1 and samples in later columns.",
    )
    parser.add_argument(
        "-g",
        "--group",
        dest="group_file",
        required=True,
        type=Path,
        help="Two-column sample/group file; comma or tab delimited.",
    )
    parser.add_argument(
        "-n",
        "--n_components",
        "--ncomponent",
        type=int,
        default=2,
        help="Number of principal components to calculate.",
    )
    parser.add_argument(
        "-l",
        "--label",
        action="store_true",
        help="Add sample IDs to the PCA plot.",
    )
    parser.add_argument(
        "-c",
        "--marker",
        choices=("o", ".", "^", "s", "D", "x"),
        default="o",
        help="Marker used in the PCA plot.",
    )
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=0.7,
        help="Point opacity between 0 and 1.",
    )
    parser.add_argument(
        "--legend_location",
        choices=(
            "best",
            "upper right",
            "lower right",
            "lower left",
            "upper left",
        ),
        default="best",
        help="Legend location.",
    )
    parser.add_argument(
        "--loading",
        action="store_true",
        help="Write the PCA loading matrix.",
    )
    parser.add_argument(
        "--no_standardize",
        action="store_true",
        help="Do not standardize CpGs before PCA.",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=8.0,
        help="Plot width in inches.",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=8.0,
        help="Plot height in inches.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Plot resolution in dots per inch.",
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        "--output",
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


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run PCA and write scores, optional loadings, and a scatter plot."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.n_components < 2:
        parser.error("--n_components must be at least 2")
    if not 0.0 <= args.alpha <= 1.0:
        parser.error("--alpha must be between 0 and 1")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be positive")
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")

    try:
        printlog(f'Reading input file: "{args.input_file}" ...')
        beta = read_beta_matrix(args.input_file)

        printlog(f'Reading group file: "{args.group_file}" ...')
        groups = read_group_file(args.group_file)

        beta_samples = list(beta.columns[1:])
        common_samples = [sample for sample in beta_samples if sample in groups.index]
        if len(common_samples) < 2:
            raise PCAError(
                "Fewer than two samples are shared by the beta matrix and group file."
            )

        missing_from_groups = [sample for sample in beta_samples if sample not in groups.index]
        if missing_from_groups:
            printlog(
                "Warning: excluding samples without group assignments: "
                + ", ".join(missing_from_groups[:10])
            )

        matrix = beta.set_index("CpG_ID")[common_samples]
        rows_before = matrix.shape[0]
        matrix = matrix.dropna(axis=0, how="any")
        removed = rows_before - matrix.shape[0]
        printlog(f"Removed {removed} CpGs containing missing values.")

        if matrix.shape[0] < 2:
            raise PCAError("Fewer than two complete CpGs remain after filtering.")

        sample_matrix = matrix.T
        maximum_components = min(sample_matrix.shape[0], sample_matrix.shape[1])
        if args.n_components > maximum_components:
            raise PCAError(
                f"--n_components={args.n_components} exceeds the maximum allowed "
                f"value ({maximum_components}) for this matrix."
            )

        values = sample_matrix.to_numpy(dtype=float)
        if not args.no_standardize:
            printlog("Standardizing CpG values before PCA ...")
            values = StandardScaler().fit_transform(values)

        printlog(
            f"Running PCA with {sample_matrix.shape[0]} samples and "
            f"{sample_matrix.shape[1]} CpGs ..."
        )
        model = PCA(n_components=args.n_components, random_state=0)
        transformed = model.fit_transform(values)

        pc_names = [f"PC{index + 1}" for index in range(args.n_components)]
        scores = pd.DataFrame(
            transformed,
            index=sample_matrix.index,
            columns=pc_names,
        )
        scores.index.name = "Sample_ID"
        scores["Group_ID"] = groups.loc[scores.index].astype(str)

        for index, ratio in enumerate(model.explained_variance_ratio_, start=1):
            printlog(f"Variance explained by PC{index}: {ratio * 100:.4f}%")

        args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
        scores_path = Path(f"{args.out_prefix}.PCA.tsv")
        plot_path = Path(f"{args.out_prefix}.PCA.png")
        variance_path = Path(f"{args.out_prefix}.PCA_variance.tsv")

        scores.to_csv(scores_path, sep="\t", float_format="%.8g")
        variance = pd.DataFrame(
            {
                "component": pc_names,
                "explained_variance_ratio": model.explained_variance_ratio_,
                "explained_variance_percent": model.explained_variance_ratio_ * 100,
                "cumulative_variance_percent": np.cumsum(
                    model.explained_variance_ratio_
                )
                * 100,
            }
        )
        variance.to_csv(variance_path, sep="\t", index=False, float_format="%.8g")

        if args.loading:
            loadings_path = Path(f"{args.out_prefix}.PCA_loadings.tsv")
            loadings = pd.DataFrame(
                model.components_.T,
                index=sample_matrix.columns,
                columns=pc_names,
            )
            loadings.index.name = "CpG_ID"
            loadings.to_csv(loadings_path, sep="\t", float_format="%.8g")

        plot_pca(
            scores=scores,
            explained_variance=model.explained_variance_ratio_,
            output_path=plot_path,
            label_samples=args.label,
            marker=args.marker,
            alpha=args.alpha,
            legend_location=args.legend_location,
            width=args.width,
            height=args.height,
            dpi=args.dpi,
        )

        printlog(f'PCA scores written to "{scores_path}".')
        printlog(f'PCA variance written to "{variance_path}".')
        printlog(f'PCA plot written to "{plot_path}".')

    except PCAError as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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

"""Perform UMAP analysis on DNA methylation samples.

The input matrix must contain CpG identifiers in the first column and sample
beta values in the remaining columns. CpGs are treated as features and samples
as observations.

An optional sample-group file can be supplied to color samples in the plot.
Group files may be comma- or tab-delimited, with or without a header.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class UMAPError(RuntimeError):
    """Raised when input validation or UMAP analysis fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_UMAP",
        description="Perform UMAP dimensionality reduction on DNA methylation samples.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input beta-value matrix. The first column contains CpG IDs and "
            "remaining columns contain samples. Compressed input is supported."
        ),
    )
    parser.add_argument(
        "-g",
        "--group",
        dest="group_file",
        type=Path,
        default=None,
        help=(
            "Optional CSV/TSV file containing sample IDs and group labels. "
            "Headered and headerless files are supported."
        ),
    )
    parser.add_argument(
        "-n",
        "--n_components",
        "--ncomponent",
        dest="n_components",
        type=int,
        choices=(2, 3),
        default=2,
        help="Number of UMAP embedding dimensions.",
    )
    parser.add_argument(
        "--n_neighbors",
        "--nneighbors",
        dest="n_neighbors",
        type=int,
        default=15,
        help=(
            "Number of neighboring samples used to learn local manifold structure. "
            "Must be at least 2 and smaller than the number of samples."
        ),
    )
    parser.add_argument(
        "--min_dist",
        "--min-dist",
        dest="min_dist",
        type=float,
        default=0.2,
        help="Minimum distance between points in the low-dimensional embedding.",
    )
    parser.add_argument(
        "--spread",
        type=float,
        default=1.0,
        help="Effective scale of embedded points.",
    )
    parser.add_argument(
        "--metric",
        default="euclidean",
        help="Distance metric accepted by umap-learn.",
    )
    parser.add_argument(
        "--init",
        choices=("spectral", "random"),
        default="spectral",
        help="Initialization method.",
    )
    parser.add_argument(
        "--densmap",
        action="store_true",
        help="Enable densMAP to better preserve local density structure.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=99,
        help="Random seed.",
    )
    parser.add_argument(
        "--na_policy",
        choices=("drop", "impute"),
        default="drop",
        help=(
            "Missing-value handling: 'drop' removes CpGs with any missing value; "
            "'impute' replaces missing values with the per-CpG mean."
        ),
    )
    parser.add_argument(
        "--no_standardize",
        action="store_true",
        help="Do not standardize CpG features across samples before UMAP.",
    )
    parser.add_argument(
        "--label",
        action="store_true",
        help="Add sample IDs to the plot.",
    )
    parser.add_argument(
        "--marker",
        choices=("circle", "dot"),
        default="circle",
        help="Plot marker style.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.8,
        help="Point opacity.",
    )
    parser.add_argument(
        "--point_size",
        type=float,
        default=45.0,
        help="Point size.",
    )
    parser.add_argument(
        "--legend_location",
        choices=("best", "upper right", "lower right", "lower left", "upper left"),
        default="best",
        help="Legend location.",
    )
    parser.add_argument(
        "--format",
        choices=("png", "pdf", "both"),
        default="pdf",
        help="Plot output format.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution of PNG output.",
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
        "--no_plot",
        action="store_true",
        help="Write UMAP coordinates without generating a plot.",
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
    if args.group_file is not None and not args.group_file.is_file():
        parser.error(f"Group file does not exist: {args.group_file}")
    if args.n_neighbors < 2:
        parser.error("--n_neighbors must be at least 2")
    if not 0.0 <= args.min_dist < 1.0:
        parser.error("--min_dist must be in the interval [0, 1)")
    if args.spread <= 0:
        parser.error("--spread must be greater than zero")
    if args.min_dist > args.spread:
        parser.error("--min_dist cannot be greater than --spread")
    if not 0.0 < args.alpha <= 1.0:
        parser.error("--alpha must be in the interval (0, 1]")
    if args.point_size <= 0:
        parser.error("--point_size must be greater than zero")
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def read_beta_matrix(path: Path, na_policy: str) -> pd.DataFrame:
    """Read and preprocess the beta-value matrix."""
    printlog(f'Reading beta-value matrix: "{path}"')

    try:
        data = pd.read_csv(
            path,
            sep=None,
            engine="python",
            compression="infer",
            index_col=0,
        )
    except Exception as exc:
        raise UMAPError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] < 3:
        raise UMAPError("At least three samples are recommended for UMAP.")

    data.index = data.index.astype(str)
    data.columns = data.columns.astype(str)

    if data.index.duplicated().any():
        duplicates = data.index[data.index.duplicated()].unique().tolist()
        raise UMAPError(
            "Duplicate CpG IDs were found: " + ", ".join(duplicates[:10])
        )
    if data.columns.duplicated().any():
        duplicates = data.columns[data.columns.duplicated()].tolist()
        raise UMAPError(
            "Duplicate sample IDs were found: " + ", ".join(duplicates)
        )

    data = data.apply(pd.to_numeric, errors="coerce")

    finite = data.to_numpy(dtype=float)
    out_of_range = np.isfinite(finite) & ((finite < 0.0) | (finite > 1.0))
    if out_of_range.any():
        printlog(
            f"Warning: {int(out_of_range.sum())} finite values fall outside "
            "the expected beta-value range [0, 1]."
        )

    if na_policy == "drop":
        before = len(data)
        data = data.dropna(axis=0, how="any")
        printlog(f"Removed {before - len(data)} CpGs containing missing values.")
    else:
        all_missing = data.isna().all(axis=1)
        if all_missing.any():
            data = data.loc[~all_missing]
            printlog(f"Removed {int(all_missing.sum())} all-missing CpGs.")

        row_means = data.mean(axis=1)
        data = data.T.fillna(row_means).T
        printlog("Imputed missing values using per-CpG means.")

    if data.empty:
        raise UMAPError("No CpGs remain after missing-value handling.")

    zero_variance = data.var(axis=1, ddof=0) == 0
    if zero_variance.any():
        data = data.loc[~zero_variance]
        printlog(f"Removed {int(zero_variance.sum())} zero-variance CpGs.")

    if data.empty:
        raise UMAPError("No variable CpGs remain for analysis.")

    printlog(f"Samples used: {data.shape[1]}")
    printlog(f"CpGs used: {data.shape[0]}")
    return data


def read_group_file(path: Optional[Path], sample_names: list[str]) -> pd.DataFrame:
    """Read, validate, and align sample-group metadata."""
    if path is None:
        return pd.DataFrame(
            {"Group": ["All_samples"] * len(sample_names)},
            index=pd.Index(sample_names, name="Sample_ID"),
        )

    printlog(f'Reading group file: "{path}"')

    try:
        group = pd.read_csv(
            path,
            sep=None,
            engine="python",
            header=None,
            dtype=str,
            compression="infer",
        )
    except Exception as exc:
        raise UMAPError(f"Cannot read group file {path}: {exc}") from exc

    if group.shape[1] < 2:
        raise UMAPError("Group file must contain at least two columns.")

    group = group.iloc[:, :2].copy()
    group.columns = ["Sample_ID", "Group"]
    group["Sample_ID"] = group["Sample_ID"].str.strip().str.lstrip("\ufeff")
    group["Group"] = group["Group"].str.strip()

    if not group.empty:
        first_sample = str(group.iloc[0, 0]).lower()
        first_group = str(group.iloc[0, 1]).lower()
        if first_sample in {"sample", "sample_id", "sampleid"} and first_group in {
            "group", "group_id", "groupid", "class", "label"
        }:
            group = group.iloc[1:].copy()

    group = group.dropna(subset=["Sample_ID", "Group"])

    if group["Sample_ID"].duplicated().any():
        duplicates = group.loc[
            group["Sample_ID"].duplicated(keep=False), "Sample_ID"
        ].unique()
        raise UMAPError(
            "Duplicate sample IDs in group file: " + ", ".join(duplicates)
        )

    group = group.set_index("Sample_ID")
    missing = [sample for sample in sample_names if sample not in group.index]
    if missing:
        raise UMAPError(
            "Samples missing from group file: " + ", ".join(missing)
        )

    return group.loc[sample_names, ["Group"]]


def run_umap(beta_data: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Run UMAP and return sample coordinates."""
    sample_data = beta_data.T
    n_samples = sample_data.shape[0]

    if args.n_neighbors >= n_samples:
        raise UMAPError(
            f"--n_neighbors ({args.n_neighbors}) must be smaller than "
            f"the number of samples ({n_samples})."
        )

    matrix = sample_data.to_numpy(dtype=float)

    if not args.no_standardize:
        printlog("Standardizing CpG features across samples ...")
        matrix = StandardScaler().fit_transform(matrix)

    printlog(
        "Running UMAP with "
        f"{n_samples} samples, {sample_data.shape[1]} CpGs, "
        f"n_neighbors={args.n_neighbors}, min_dist={args.min_dist} ..."
    )

    reducer = umap.UMAP(
        n_components=args.n_components,
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        spread=args.spread,
        metric=args.metric,
        init=args.init,
        densmap=args.densmap,
        random_state=args.seed,
    )

    embedding = reducer.fit_transform(matrix)
    columns = [f"UMAP{i + 1}" for i in range(args.n_components)]
    return pd.DataFrame(embedding, index=sample_data.index, columns=columns)


def plot_embedding(
    results: pd.DataFrame,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """Create a two-dimensional UMAP scatter plot."""
    groups = results["Group"].astype(str)
    unique_groups = list(dict.fromkeys(groups))
    marker = "o" if args.marker == "circle" else "."

    figure, axis = plt.subplots(figsize=(args.width, args.height))

    for group_name in unique_groups:
        mask = groups == group_name
        axis.scatter(
            results.loc[mask, "UMAP1"],
            results.loc[mask, "UMAP2"],
            label=group_name,
            marker=marker,
            alpha=args.alpha,
            s=args.point_size,
        )

    if args.label:
        for sample_name, row in results.iterrows():
            axis.annotate(
                str(sample_name),
                (row["UMAP1"], row["UMAP2"]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )

    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title("UMAP of DNA methylation samples")

    if len(unique_groups) > 1 or unique_groups != ["All_samples"]:
        axis.legend(loc=args.legend_location, title="Group")

    figure.tight_layout()
    save_kwargs = {"dpi": args.dpi} if output_path.suffix.lower() == ".png" else {}
    figure.savefig(output_path, **save_kwargs)
    plt.close(figure)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        beta_data = read_beta_matrix(args.input_file, args.na_policy)
        coordinates = run_umap(beta_data, args)
        metadata = read_group_file(args.group_file, coordinates.index.tolist())

        results = coordinates.join(metadata)
        results.index.name = "Sample_ID"

        result_path = Path(f"{args.out_prefix}.UMAP.tsv")
        printlog(f'Writing UMAP coordinates: "{result_path}"')
        results.to_csv(result_path, sep="\t", index=True, float_format="%.6f")

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)
            for extension in extensions:
                plot_path = Path(f"{args.out_prefix}.UMAP.{extension}")
                printlog(f'Writing UMAP plot: "{plot_path}"')
                plot_embedding(results, plot_path, args)

    except UMAPError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

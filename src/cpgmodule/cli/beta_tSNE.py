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

"""Perform t-SNE analysis on DNA methylation samples.

The input matrix must contain CpG identifiers in the first column and sample
beta values in the remaining columns. CpGs are treated as features and samples
as observations.

The program can optionally use a sample-group file to color samples in the
plot. Group files may be comma- or tab-delimited and may contain either:

    Sample,Group
    Sample_01,normal
    Sample_02,tumor

or headerless rows:

    Sample_01,normal
    Sample_02,tumor
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class TSNEError(RuntimeError):
    """Raised when input validation or t-SNE analysis fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_tSNE",
        description="Perform t-SNE analysis on DNA methylation samples.",
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
            "Optional CSV or TSV file containing sample IDs and group labels. "
            "Headered and headerless files are supported."
        ),
    )
    parser.add_argument(
        "-p",
        "--perplexity",
        type=float,
        default=5.0,
        help=(
            "t-SNE perplexity. It must be greater than zero and smaller than "
            "the number of samples."
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
        help="Number of t-SNE embedding dimensions.",
    )
    parser.add_argument(
        "--max_iter",
        "--n_iter",
        dest="max_iter",
        type=int,
        default=5000,
        help="Maximum number of optimization iterations.",
    )
    parser.add_argument(
        "--learning_rate",
        default="auto",
        help=(
            "Learning rate. Use 'auto' or a positive numeric value."
        ),
    )
    parser.add_argument(
        "--init",
        choices=("pca", "random"),
        default="pca",
        help="Initialization method.",
    )
    parser.add_argument(
        "--metric",
        default="euclidean",
        help="Distance metric accepted by sklearn.manifold.TSNE.",
    )
    parser.add_argument(
        "--early_exaggeration",
        type=float,
        default=12.0,
        help="Early-exaggeration factor.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed.",
    )
    parser.add_argument(
        "--no_standardize",
        action="store_true",
        help="Do not standardize CpG features before t-SNE.",
    )
    parser.add_argument(
        "--na_policy",
        choices=("drop", "impute"),
        default="drop",
        help=(
            "'drop' removes CpGs containing any missing value; 'impute' "
            "replaces missing values with the per-CpG mean."
        ),
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
        help="Write coordinates without generating a plot.",
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
    if args.perplexity <= 0:
        parser.error("--perplexity must be greater than zero")
    if args.max_iter < 250:
        parser.error("--max_iter must be at least 250")
    if args.early_exaggeration <= 0:
        parser.error("--early_exaggeration must be greater than zero")
    if not 0.0 < args.alpha <= 1.0:
        parser.error("--alpha must be in the interval (0, 1]")
    if args.point_size <= 0:
        parser.error("--point_size must be greater than zero")
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")

    if args.learning_rate != "auto":
        try:
            learning_rate = float(args.learning_rate)
        except ValueError:
            parser.error("--learning_rate must be 'auto' or a positive number")
        if learning_rate <= 0:
            parser.error("--learning_rate must be greater than zero")
        args.learning_rate = learning_rate


def read_beta_matrix(path: Path, na_policy: str) -> pd.DataFrame:
    """Read and preprocess a beta-value matrix."""
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
        raise TSNEError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] < 2:
        raise TSNEError("At least two samples are required for t-SNE.")

    data.index = data.index.astype(str)
    data.columns = data.columns.astype(str)

    if data.index.duplicated().any():
        duplicates = data.index[data.index.duplicated()].unique().tolist()
        raise TSNEError(
            "Duplicate CpG IDs were found: " + ", ".join(duplicates[:10])
        )

    if data.columns.duplicated().any():
        duplicates = data.columns[data.columns.duplicated()].tolist()
        raise TSNEError(
            "Duplicate sample IDs were found: " + ", ".join(duplicates)
        )

    data = data.apply(pd.to_numeric, errors="coerce")

    if na_policy == "drop":
        before = len(data)
        data = data.dropna(axis=0, how="any")
        printlog(f"Removed {before - len(data)} CpGs containing missing values.")
    else:
        all_missing = data.isna().all(axis=1)
        if all_missing.any():
            data = data.loc[~all_missing]
            printlog(f"Removed {int(all_missing.sum())} all-missing CpGs.")
        data = data.T.fillna(data.mean(axis=1)).T
        printlog("Imputed missing values using per-CpG means.")

    if data.empty:
        raise TSNEError("No CpGs remain after missing-value handling.")

    zero_variance = data.var(axis=1, ddof=0) == 0
    if zero_variance.any():
        data = data.loc[~zero_variance]
        printlog(f"Removed {int(zero_variance.sum())} zero-variance CpGs.")

    if data.empty:
        raise TSNEError("No variable CpGs remain for analysis.")

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
        raise TSNEError(f"Cannot read group file {path}: {exc}") from exc

    if group.shape[1] < 2:
        raise TSNEError("Group file must contain at least two columns.")

    group = group.iloc[:, :2].copy()
    group.columns = ["Sample_ID", "Group"]
    group["Sample_ID"] = group["Sample_ID"].str.strip().str.lstrip("\ufeff")
    group["Group"] = group["Group"].str.strip()

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

    group = group.dropna(subset=["Sample_ID", "Group"])

    if group["Sample_ID"].duplicated().any():
        duplicates = group.loc[
            group["Sample_ID"].duplicated(keep=False), "Sample_ID"
        ].unique()
        raise TSNEError(
            "Duplicate sample IDs in group file: " + ", ".join(duplicates)
        )

    group = group.set_index("Sample_ID")
    missing = [sample for sample in sample_names if sample not in group.index]
    if missing:
        raise TSNEError(
            "Samples missing from group file: " + ", ".join(missing)
        )

    return group.loc[sample_names, ["Group"]]


def run_tsne(
    beta_data: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Run t-SNE and return sample coordinates."""
    sample_data = beta_data.T
    n_samples = sample_data.shape[0]

    if args.perplexity >= n_samples:
        raise TSNEError(
            f"--perplexity ({args.perplexity}) must be smaller than the "
            f"number of samples ({n_samples})."
        )

    matrix = sample_data.to_numpy(dtype=float)

    if not args.no_standardize:
        printlog("Standardizing CpG features across samples ...")
        matrix = StandardScaler().fit_transform(matrix)

    printlog(
        "Running t-SNE with "
        f"{n_samples} samples, {sample_data.shape[1]} CpGs, "
        f"perplexity={args.perplexity} ..."
    )

    model = TSNE(
        n_components=args.n_components,
        perplexity=args.perplexity,
        learning_rate=args.learning_rate,
        max_iter=args.max_iter,
        init=args.init,
        metric=args.metric,
        early_exaggeration=args.early_exaggeration,
        random_state=args.seed,
    )

    embedding = model.fit_transform(matrix)
    coordinate_names = [f"tSNE{i + 1}" for i in range(args.n_components)]

    return pd.DataFrame(
        embedding,
        index=sample_data.index,
        columns=coordinate_names,
    )


def plot_embedding(
    results: pd.DataFrame,
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """Create a two-dimensional t-SNE scatter plot."""
    if args.n_components < 2:
        raise TSNEError("At least two dimensions are required for plotting.")

    groups = results["Group"].astype(str)
    unique_groups = list(dict.fromkeys(groups))
    marker = "o" if args.marker == "circle" else "."

    figure, axis = plt.subplots(figsize=(args.width, args.height))

    for group_name in unique_groups:
        mask = groups == group_name
        axis.scatter(
            results.loc[mask, "tSNE1"],
            results.loc[mask, "tSNE2"],
            label=group_name,
            marker=marker,
            alpha=args.alpha,
            s=args.point_size,
        )

    if args.label:
        for sample_name, row in results.iterrows():
            axis.annotate(
                str(sample_name),
                (row["tSNE1"], row["tSNE2"]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=7,
            )

    axis.set_xlabel("tSNE1")
    axis.set_ylabel("tSNE2")
    axis.set_title("t-SNE of DNA methylation samples")

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
        coordinates = run_tsne(beta_data, args)
        metadata = read_group_file(args.group_file, coordinates.index.tolist())

        results = coordinates.join(metadata)
        results.index.name = "Sample_ID"

        result_path = Path(f"{args.out_prefix}.tSNE.tsv")
        printlog(f'Writing t-SNE coordinates: "{result_path}"')
        results.to_csv(
            result_path,
            sep="\t",
            index=True,
            float_format="%.6f",
        )

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)
            for extension in extensions:
                plot_path = Path(f"{args.out_prefix}.tSNE.{extension}")
                printlog(f'Writing t-SNE plot: "{plot_path}"')
                plot_embedding(results, plot_path, args)

    except TSNEError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

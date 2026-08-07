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

"""Create stacked bar plots of DNA methylation beta-value distributions.

For each sample, beta values are grouped into four intervals:

* [0.00, 0.25]
* (0.25, 0.50]
* (0.50, 0.75]
* (0.75, 1.00]

The program writes a table containing per-sample counts and percentages, and
optionally generates a stacked bar plot.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


BIN_LABELS = (
    "beta [0.00 - 0.25]",
    "beta (0.25 - 0.50]",
    "beta (0.50 - 0.75]",
    "beta (0.75 - 1.00]",
)

BIN_EDGES = (-np.inf, 0.25, 0.50, 0.75, np.inf)


class StackedBarplotError(RuntimeError):
    """Raised when input validation or plotting fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_stacked_barplot",
        description=(
            "Summarize beta-value distributions for each sample and create "
            "a stacked bar plot."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input beta-value matrix. The first column contains CpG IDs and "
            "remaining columns contain sample beta values. Compressed input "
            "is supported."
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
        default=10.0,
        help="Plot width in inches.",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=6.0,
        help="Plot height in inches.",
    )
    parser.add_argument(
        "--no_plot",
        action="store_true",
        help="Write the summary table without generating a plot.",
    )
    parser.add_argument(
        "--allow_out_of_range",
        action="store_true",
        help=(
            "Ignore beta values outside [0, 1]. By default, the program exits "
            "if finite out-of-range values are detected."
        ),
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
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def read_beta_matrix(path: Path) -> pd.DataFrame:
    """Read and validate a beta-value matrix."""
    printlog(f'Reading beta-value matrix: "{path}"')

    try:
        data = pd.read_csv(
            path,
            sep=None,
            engine="python",
            compression="infer",
        )
    except Exception as exc:
        raise StackedBarplotError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] < 2:
        raise StackedBarplotError(
            "Input must contain one CpG-ID column and at least one sample column."
        )

    if data.columns[1:].duplicated().any():
        duplicates = data.columns[1:][data.columns[1:].duplicated()].tolist()
        raise StackedBarplotError(
            "Duplicate sample IDs were found: " + ", ".join(map(str, duplicates))
        )

    return data


def summarize_beta_values(
    data: pd.DataFrame,
    allow_out_of_range: bool = False,
) -> pd.DataFrame:
    """Calculate beta-value counts and percentages for each sample."""
    sample_names = list(data.columns[1:])
    summaries = []

    for sample in sample_names:
        values = pd.to_numeric(data[sample], errors="coerce")
        finite_values = values[np.isfinite(values)]

        out_of_range = finite_values[(finite_values < 0.0) | (finite_values > 1.0)]
        if not out_of_range.empty and not allow_out_of_range:
            raise StackedBarplotError(
                f'Sample "{sample}" contains {len(out_of_range)} finite beta '
                "values outside [0, 1]. Use --allow_out_of_range to ignore them."
            )

        valid = finite_values[(finite_values >= 0.0) & (finite_values <= 1.0)]
        missing_count = int(len(values) - len(finite_values))
        out_of_range_count = int(len(out_of_range))

        if valid.empty:
            raise StackedBarplotError(
                f'Sample "{sample}" contains no valid beta values in [0, 1].'
            )

        categories = pd.cut(
            valid,
            bins=BIN_EDGES,
            labels=BIN_LABELS,
            include_lowest=True,
            right=True,
        )
        counts = categories.value_counts(sort=False).astype(int)
        total_valid = int(counts.sum())

        for label in BIN_LABELS:
            count = int(counts.get(label, 0))
            summaries.append(
                {
                    "sample": sample,
                    "beta_range": label,
                    "count": count,
                    "percentage": 100.0 * count / total_valid,
                    "valid_values": total_valid,
                    "missing_values": missing_count,
                    "out_of_range_values": out_of_range_count,
                }
            )

    return pd.DataFrame(summaries)


def plot_summary(
    summary: pd.DataFrame,
    output_path: Path,
    width: float,
    height: float,
    dpi: int,
) -> None:
    """Generate a stacked bar plot from the summary table."""
    percentage_table = summary.pivot(
        index="sample",
        columns="beta_range",
        values="percentage",
    ).reindex(columns=BIN_LABELS)

    sample_names = percentage_table.index.tolist()
    x = np.arange(len(sample_names))
    bottom = np.zeros(len(sample_names), dtype=float)

    figure, axis = plt.subplots(figsize=(width, height))

    for label in BIN_LABELS:
        values = percentage_table[label].to_numpy(dtype=float)
        axis.bar(
            x,
            values,
            bottom=bottom,
            width=0.85,
            edgecolor="white",
            linewidth=0.4,
            label=label,
        )
        bottom += values

    axis.set_ylabel("Percentage of valid CpGs")
    axis.set_xlabel("Sample")
    axis.set_ylim(0.0, 100.0)
    axis.set_xticks(x)
    axis.set_xticklabels(sample_names, rotation=90)
    axis.set_title("DNA methylation beta-value distribution")
    axis.legend(
        title="Beta range",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
    )

    figure.tight_layout()
    save_kwargs = {"dpi": dpi} if output_path.suffix.lower() == ".png" else {}
    figure.savefig(output_path, **save_kwargs)
    plt.close(figure)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        beta_data = read_beta_matrix(args.input_file)
        summary = summarize_beta_values(
            beta_data,
            allow_out_of_range=args.allow_out_of_range,
        )

        summary_path = Path(f"{args.out_prefix}.stacked_barplot.tsv")
        printlog(f'Writing beta-distribution summary: "{summary_path}"')
        summary.to_csv(
            summary_path,
            sep="\t",
            index=False,
            float_format="%.6f",
        )

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)
            for extension in extensions:
                plot_path = Path(f"{args.out_prefix}.stacked_barplot.{extension}")
                printlog(f'Writing stacked bar plot: "{plot_path}"')
                plot_summary(
                    summary=summary,
                    output_path=plot_path,
                    width=args.width,
                    height=args.height,
                    dpi=args.dpi,
                )

    except StackedBarplotError as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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

"""Summarize and visualize CpG distributions across chromosomes.

The program accepts one or more BED3+ CpG files and a chromosome-size file.
For each sample, it reports:

* total CpG count per chromosome
* percentage of all CpGs assigned to each chromosome
* CpG density per megabase

Chromosomes and their plotting order are determined by the chromosome-size file.
Compressed CpG input files are supported through the CpGtools reader.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class ChromDistributionError(RuntimeError):
    """Raised when chromosome distribution analysis fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_distrb_chrom",
        description="Calculate CpG distributions across chromosomes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_files",
        required=True,
        nargs="+",
        type=Path,
        help=(
            "One or more CpG BED3+ files. Files may also be supplied as a "
            "single comma-separated argument for backward compatibility."
        ),
    )
    parser.add_argument(
        "-n",
        "--names",
        nargs="+",
        default=None,
        help=(
            "Optional sample names corresponding to input files. Names may also "
            "be supplied as one comma-separated argument."
        ),
    )
    parser.add_argument(
        "-s",
        "--chrom_size",
        "--chrom-size",
        dest="chrom_size",
        required=True,
        type=Path,
        help=(
            "Chromosome-size file with chromosome name in column 1 and size "
            "in column 2. Chromosome order in this file determines output order."
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
        default=12.0,
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
        help="Write summary tables without generating plots.",
    )
    parser.add_argument(
        "--max_plot_samples",
        type=int,
        default=12,
        help=(
            "Maximum number of samples to include in grouped bar plots. "
            "Set to 0 to disable this limit."
        ),
    )
    parser.add_argument(
        "--strict_chromosomes",
        action="store_true",
        help=(
            "Exit if CpG files contain chromosomes not present in the "
            "chromosome-size file. By default, such chromosomes are ignored "
            "with a warning."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def expand_comma_separated(values: Sequence[object]) -> list[str]:
    """Expand list arguments that may contain comma-separated values."""
    expanded: list[str] = []
    for value in values:
        expanded.extend(
            item.strip()
            for item in str(value).split(",")
            if item.strip()
        )
    return expanded


def validate_args(args: argparse.Namespace, parser: argparse.ArgumentParser) -> None:
    """Normalize and validate command-line arguments."""
    args.input_files = [Path(p) for p in expand_comma_separated(args.input_files)]

    if not args.input_files:
        parser.error("At least one input CpG file is required.")

    for path in args.input_files:
        if not path.is_file():
            parser.error(f"Input CpG file does not exist: {path}")

    if not args.chrom_size.is_file():
        parser.error(f"Chromosome-size file does not exist: {args.chrom_size}")

    if args.names is not None:
        args.names = expand_comma_separated(args.names)
        if len(args.names) != len(args.input_files):
            parser.error(
                "--names must contain exactly one name for each input CpG file."
            )

    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")

    if args.dpi < 1:
        parser.error("--dpi must be at least 1")

    if args.max_plot_samples < 0:
        parser.error("--max_plot_samples must be non-negative")


def read_chrom_sizes(path: Path) -> pd.DataFrame:
    """Read chromosome names and sizes while preserving file order."""
    printlog(f'Reading chromosome-size file: "{path}"')

    records: list[tuple[str, int]] = []
    seen: set[str] = set()

    for line_number, raw_line in enumerate(ireader.reader(str(path)), start=1):
        line = raw_line.strip()

        if not line or line.startswith(("#", "track", "browser")):
            continue

        fields = line.split()
        if len(fields) < 2:
            raise ChromDistributionError(
                f"Chromosome-size line {line_number} has fewer than two columns."
            )

        chrom = fields[0]

        try:
            size = int(fields[1])
        except ValueError as exc:
            raise ChromDistributionError(
                f"Invalid chromosome size on line {line_number}: {line}"
            ) from exc

        if size <= 0:
            raise ChromDistributionError(
                f"Chromosome size must be positive on line {line_number}: {line}"
            )

        if chrom in seen:
            raise ChromDistributionError(
                f"Duplicate chromosome {chrom!r} in chromosome-size file."
            )

        seen.add(chrom)
        records.append((chrom, size))

    if not records:
        raise ChromDistributionError(
            "Chromosome-size file contains no valid chromosome records."
        )

    return pd.DataFrame(records, columns=["Chromosome", "ChromSize"])


def count_cpgs_by_chromosome(path: Path) -> tuple[dict[str, int], int]:
    """Count valid BED records by chromosome."""
    printlog(f'Reading CpG BED file: "{path}"')

    counts: dict[str, int] = {}
    invalid = 0

    for line_number, raw_line in enumerate(ireader.reader(str(path)), start=1):
        line = raw_line.strip()

        if not line or line.startswith(("#", "track", "browser")):
            continue

        fields = line.split()
        if len(fields) < 3:
            invalid += 1
            continue

        chrom = fields[0]

        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError:
            invalid += 1
            continue

        if start < 0 or end <= start:
            invalid += 1
            continue

        counts[chrom] = counts.get(chrom, 0) + 1

    return counts, invalid


def build_summary(
    input_files: Sequence[Path],
    sample_names: Sequence[str],
    chrom_sizes: pd.DataFrame,
    strict_chromosomes: bool,
) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    """Build chromosome-level count, percentage, and density summaries."""
    summary = chrom_sizes.copy()
    allowed_chromosomes = set(summary["Chromosome"])
    ignored: dict[str, list[str]] = {}

    for path, sample in zip(input_files, sample_names):
        counts, invalid = count_cpgs_by_chromosome(path)

        if invalid:
            printlog(
                f'Warning: ignored {invalid} invalid BED records in "{path}".'
            )

        extra_chromosomes = sorted(set(counts) - allowed_chromosomes)
        if extra_chromosomes:
            ignored[sample] = extra_chromosomes
            message = (
                f'Sample "{sample}" contains CpGs on chromosomes absent from '
                f'the chromosome-size file: {", ".join(extra_chromosomes)}'
            )
            if strict_chromosomes:
                raise ChromDistributionError(message)
            printlog("Warning: " + message + ". These CpGs will be ignored.")

        count_values = np.asarray(
            [counts.get(chrom, 0) for chrom in summary["Chromosome"]],
            dtype=int,
        )
        total_count = int(count_values.sum())

        if total_count > 0:
            percentages = count_values * 100.0 / total_count
        else:
            percentages = np.zeros(len(count_values), dtype=float)

        density = (
            count_values
            * 1_000_000.0
            / summary["ChromSize"].to_numpy(dtype=float)
        )

        summary[f"{sample}.CpG_count"] = count_values
        summary[f"{sample}.CpG_percent"] = percentages
        summary[f"{sample}.CpG_perMb"] = density

    return summary, ignored


def write_summary_tables(
    summary: pd.DataFrame,
    sample_names: Sequence[str],
    out_prefix: Path,
) -> None:
    """Write combined and metric-specific summary tables."""
    combined_path = Path(f"{out_prefix}.chrom_distribution.tsv")
    printlog(f'Writing chromosome distribution table: "{combined_path}"')
    summary.to_csv(
        combined_path,
        sep="\t",
        index=False,
        float_format="%.6f",
    )

    for suffix, metric in (
        ("counts", "CpG_count"),
        ("percent", "CpG_percent"),
        ("perMb", "CpG_perMb"),
    ):
        columns = ["Chromosome", "ChromSize"] + [
            f"{sample}.{metric}" for sample in sample_names
        ]
        path = Path(f"{out_prefix}.chrom_{suffix}.tsv")
        printlog(f'Writing chromosome {suffix} table: "{path}"')
        summary.loc[:, columns].to_csv(
            path,
            sep="\t",
            index=False,
            float_format="%.6f",
        )


def plot_metric(
    summary: pd.DataFrame,
    sample_names: Sequence[str],
    metric: str,
    ylabel: str,
    title: str,
    output_path: Path,
    *,
    width: float,
    height: float,
    dpi: int,
) -> None:
    """Create a grouped chromosome bar plot for one metric."""
    chromosomes = summary["Chromosome"].astype(str).tolist()
    x = np.arange(len(chromosomes))
    n_samples = len(sample_names)

    total_width = 0.82
    bar_width = total_width / max(n_samples, 1)

    figure, axis = plt.subplots(figsize=(width, height))

    for index, sample in enumerate(sample_names):
        values = summary[f"{sample}.{metric}"].to_numpy(dtype=float)
        offset = (index - (n_samples - 1) / 2.0) * bar_width

        axis.bar(
            x + offset,
            values,
            width=bar_width,
            label=sample,
        )

    axis.set_xlabel("Chromosome")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    axis.set_xticks(x)
    axis.set_xticklabels(chromosomes, rotation=45, ha="right")

    if n_samples > 1:
        axis.legend(
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
        )

    figure.tight_layout()

    save_kwargs = {"dpi": dpi} if output_path.suffix.lower() == ".png" else {}
    figure.savefig(output_path, **save_kwargs)
    plt.close(figure)


def generate_plots(
    summary: pd.DataFrame,
    sample_names: Sequence[str],
    args: argparse.Namespace,
) -> None:
    """Generate count, percentage, and density chromosome plots."""
    extensions = ("png", "pdf") if args.format == "both" else (args.format,)

    metrics = (
        (
            "CpG_count",
            "CpG count",
            "CpG count by chromosome",
            "CpG_total",
        ),
        (
            "CpG_percent",
            "CpG percentage (%)",
            "Percentage of CpGs by chromosome",
            "CpG_percent",
        ),
        (
            "CpG_perMb",
            "CpGs per Mb",
            "CpG density by chromosome",
            "CpG_perMb",
        ),
    )

    for metric, ylabel, title, suffix in metrics:
        for extension in extensions:
            output_path = Path(
                f"{args.out_prefix}.{suffix}.{extension}"
            )
            printlog(f'Writing chromosome plot: "{output_path}"')
            plot_metric(
                summary=summary,
                sample_names=sample_names,
                metric=metric,
                ylabel=ylabel,
                title=title,
                output_path=output_path,
                width=args.width,
                height=args.height,
                dpi=args.dpi,
            )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    sample_names = (
        args.names
        if args.names is not None
        else [path.name for path in args.input_files]
    )

    if len(set(sample_names)) != len(sample_names):
        parser.error("Sample names must be unique.")

    try:
        chrom_sizes = read_chrom_sizes(args.chrom_size)

        summary, _ = build_summary(
            input_files=args.input_files,
            sample_names=sample_names,
            chrom_sizes=chrom_sizes,
            strict_chromosomes=args.strict_chromosomes,
        )

        write_summary_tables(
            summary=summary,
            sample_names=sample_names,
            out_prefix=args.out_prefix,
        )

        if not args.no_plot:
            if (
                args.max_plot_samples > 0
                and len(sample_names) > args.max_plot_samples
            ):
                printlog(
                    f"Warning: {len(sample_names)} samples exceed "
                    f"--max_plot_samples={args.max_plot_samples}; "
                    "plots will not be generated."
                )
            else:
                generate_plots(
                    summary=summary,
                    sample_names=sample_names,
                    args=args,
                )

    except ChromDistributionError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

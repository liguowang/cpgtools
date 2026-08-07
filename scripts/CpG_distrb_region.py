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

"""Summarize CpG distribution across user-defined genomic region classes.

The program accepts one CpG BED3+ file and multiple BED3+ annotation files.
The order of annotation files defines priority. Regions from lower-priority
files are merged and then have all higher-priority regions subtracted, making
the final region classes mutually exclusive.

For each annotation class, the program reports:

* priority order
* annotation name
* number of non-overlapping regions
* total annotated size in bp
* raw CpG count
* CpG density per kilobase

Compressed BED files are supported through the CpGtools reader.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd

from cpgmodule import BED
from cpgmodule._version import __version__
from cpgmodule.utils import count_over_range, printlog, read_CpG_bed, read_bed_as_list


class RegionDistributionError(RuntimeError):
    """Raised when CpG distribution analysis across regions fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_distrb_region",
        description=(
            "Calculate CpG distribution across prioritized user-defined "
            "genomic region classes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--cpg",
        dest="cpg_file",
        required=True,
        type=Path,
        help=(
            "Input CpG BED3+ file. Only chromosome/start/end are required. "
            "Compressed input is supported."
        ),
    )
    parser.add_argument(
        "-b",
        "--bed",
        dest="bed_files",
        required=True,
        nargs="+",
        type=Path,
        help=(
            "One or more BED3+ files defining genomic region classes, listed "
            "from highest to lowest priority. A comma-separated list is also "
            "accepted for backward compatibility."
        ),
    )
    parser.add_argument(
        "-n",
        "--names",
        nargs="+",
        default=None,
        help=(
            "Optional names for the region classes. Must match the number of "
            "BED files. Comma-separated names are also accepted."
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
        default=8.0,
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
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def expand_comma_separated(values: Sequence[object]) -> list[str]:
    """Expand argparse values that may contain comma-separated entries."""
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
    args.bed_files = [Path(p) for p in expand_comma_separated(args.bed_files)]

    if not args.cpg_file.is_file():
        parser.error(f"CpG BED file does not exist: {args.cpg_file}")

    if not args.bed_files:
        parser.error("At least one region BED file is required.")

    for bed_file in args.bed_files:
        if not bed_file.is_file():
            parser.error(f"Region BED file does not exist: {bed_file}")

    if args.names is not None:
        args.names = expand_comma_separated(args.names)
        if len(args.names) != len(args.bed_files):
            parser.error(
                "--names must contain exactly one name for each BED file."
            )

    if args.dpi < 1:
        parser.error("--dpi must be at least 1")

    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def safe_union(intervals):
    """Merge BED3 intervals and safely handle empty inputs."""
    if not intervals:
        return []
    return BED.unionBed3(intervals)


def safe_subtract(intervals, higher_priority):
    """Subtract higher-priority intervals while safely handling empties."""
    if not intervals:
        return []
    if not higher_priority:
        return intervals
    return BED.subtractBed3(intervals, higher_priority)


def load_and_prioritize_regions(
    bed_files: Sequence[Path],
    region_names: Sequence[str],
) -> list[tuple[str, list]]:
    """Load, merge, and prioritize region BED files."""
    prioritized: list[tuple[str, list]] = []
    higher_priority_union = []

    for priority, (bed_file, region_name) in enumerate(
        zip(bed_files, region_names)
    ):
        printlog(
            f'Reading region BED file #{priority}: "{bed_file}" '
            f'({region_name})'
        )

        intervals = read_bed_as_list(str(bed_file))
        merged = safe_union(intervals)

        if priority > 0:
            printlog(
                f'Subtracting higher-priority regions from "{region_name}" ...'
            )
            merged = safe_subtract(merged, higher_priority_union)

        prioritized.append((region_name, merged))

        if merged:
            higher_priority_union = safe_union(
                higher_priority_union + merged
            )

    return prioritized


def summarize_regions(
    prioritized_regions: Sequence[tuple[str, list]],
    cpg_ranges,
) -> pd.DataFrame:
    """Calculate region coverage and CpG density for each class."""
    rows = []

    for priority, (region_name, intervals) in enumerate(prioritized_regions):
        printlog(f'Counting CpGs in "{region_name}" ...')

        if intervals:
            size, count = count_over_range(intervals, cpg_ranges)
            size = int(size)
            count = int(count)
        else:
            size = 0
            count = 0

        density = count * 1000.0 / size if size > 0 else 0.0

        rows.append(
            {
                "Priority_order": priority,
                "Name": region_name,
                "Number_of_regions": len(intervals),
                "Size_of_regions_bp": size,
                "CpG_raw_count": count,
                "CpG_count_per_KB": density,
            }
        )

    return pd.DataFrame(rows)


def plot_summary(
    summary: pd.DataFrame,
    output_path: Path,
    *,
    dpi: int,
    width: float,
    height: float,
) -> None:
    """Create a CpG-density bar plot across prioritized region classes."""
    figure, axis = plt.subplots(figsize=(width, height))

    x = range(len(summary))
    axis.bar(
        x,
        summary["CpG_count_per_KB"].to_numpy(dtype=float),
        width=0.75,
    )

    axis.set_xticks(list(x))
    axis.set_xticklabels(
        summary["Name"].astype(str).tolist(),
        rotation=30,
        ha="right",
    )
    axis.set_ylabel("CpG count per kb")
    axis.set_xlabel("Genomic region class")
    axis.set_title("CpG density across prioritized genomic regions")

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

    region_names = (
        args.names
        if args.names is not None
        else [path.name for path in args.bed_files]
    )

    if len(set(region_names)) != len(region_names):
        parser.error("Region names must be unique.")

    table_path = Path(f"{args.out_prefix}.region_distribution.tsv")

    try:
        printlog(f'Reading CpG file: "{args.cpg_file}"')
        cpg_ranges = read_CpG_bed(str(args.cpg_file))

        prioritized_regions = load_and_prioritize_regions(
            bed_files=args.bed_files,
            region_names=region_names,
        )

        summary = summarize_regions(
            prioritized_regions=prioritized_regions,
            cpg_ranges=cpg_ranges,
        )

        printlog(f'Writing region distribution table: "{table_path}"')
        summary.to_csv(
            table_path,
            sep="\t",
            index=False,
            float_format="%.6f",
        )

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)

            for extension in extensions:
                plot_path = Path(
                    f"{args.out_prefix}.region_distribution.{extension}"
                )
                printlog(f'Writing region distribution plot: "{plot_path}"')
                plot_summary(
                    summary=summary,
                    output_path=plot_path,
                    dpi=args.dpi,
                    width=args.width,
                    height=args.height,
                )

    except RegionDistributionError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

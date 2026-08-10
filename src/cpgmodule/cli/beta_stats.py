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

"""Summarize DNA methylation values within user-defined genomic regions.

The input methylation file must be BED6 or BED6+ with beta values in column 5.
The region file must be BED3 or BED3+. All original region columns are preserved,
and the following summary columns are appended:

* CpG_count
* Minimum_beta
* Maximum_beta
* Mean_beta
* Median_beta
* Standard_deviation

Coordinates are interpreted as zero-based, half-open BED intervals.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional, Sequence

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog, read_CpG_bed, stats_over_range


STAT_COLUMNS = (
    "CpG_count",
    "Minimum_beta",
    "Maximum_beta",
    "Mean_beta",
    "Median_beta",
    "Standard_deviation",
)

def build_header(region_file: Path) -> list[str]:
    """Infer the number of input columns and append statistic columns."""
    for _, line in iter_region_lines(region_file):
        column_count = len(line.split())

        if column_count < 3:
            raise BetaStatsError(
                "Region file must contain at least three columns."
            )

        input_columns = ["chrom", "start", "end"]

        if column_count > 3:
            input_columns.extend(
                f"column_{index}"
                for index in range(4, column_count + 1)
            )

        return input_columns + list(STAT_COLUMNS)

    raise BetaStatsError("Region file contains no valid data records.")
    
class BetaStatsError(RuntimeError):
    """Raised when input validation or regional summarization fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_stats",
        description=(
            "Calculate CpG counts and methylation summary statistics within "
            "user-defined genomic regions. Six additional columns will be \
            appended to the input file, including: CpG_count, Minimum_beta, \
            Maximum_beta, Mean_beta, Median_beta, Standard_deviation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input BED6+ methylation file. Column 5 must contain beta values. "
            "Compressed input is supported by the CpGtools reader."
        ),
    )
    parser.add_argument(
        "-r",
        "--region",
        required=True,
        type=Path,
        help=(
            "BED3 file containing genomic regions. Original columns are "
            "preserved in the output."
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
        "--header",
        action="store_true",
        help="Write a header row to the output table.",
    )
    parser.add_argument(
        "--na_rep",
        default="NA",
        help="Text used for unavailable statistics.",
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
        parser.error(f"Input methylation file does not exist: {args.input_file}")
    if not args.region.is_file():
        parser.error(f"Region file does not exist: {args.region}")
    if "\t" in args.na_rep or "\n" in args.na_rep:
        parser.error("--na_rep cannot contain tabs or newlines")


def format_stat(value, na_rep: str) -> str:
    """Format a statistic for tabular output."""
    if value is None:
        return na_rep

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)

    if math.isnan(numeric) or math.isinf(numeric):
        return na_rep

    if numeric.is_integer():
        return str(int(numeric))

    return f"{numeric:.6f}"


def iter_region_lines(path: Path):
    """Yield nonempty, non-directive lines from a BED-like region file."""
    for line_number, raw_line in enumerate(ireader.reader(str(path)), start=1):
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()

        if not stripped:
            continue
        if stripped.startswith(("#", "track", "browser")):
            continue

        yield line_number, line


def summarize_regions(
    cpg_ranges,
    region_file: Path,
    output_file: Path,
    write_header: bool,
    na_rep: str,
) -> tuple[int, int]:
    """Summarize methylation values for each region and write the result."""
    processed = 0
    invalid = 0

    with output_file.open("w", encoding="utf-8") as output_handle:
        if write_header:
            header = build_header(region_file)
            output_handle.write("\t".join(header) + "\n")

        for line_number, line in iter_region_lines(region_file):
            fields = line.split()

            if len(fields) < 3:
                invalid += 1
                output_handle.write(
                    line + "\t" + "\t".join([na_rep] * len(STAT_COLUMNS)) + "\n"
                )
                continue

            chrom = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                invalid += 1
                output_handle.write(
                    line + "\t" + "\t".join([na_rep] * len(STAT_COLUMNS)) + "\n"
                )
                continue

            if start < 0 or end <= start:
                invalid += 1
                output_handle.write(
                    line + "\t" + "\t".join([na_rep] * len(STAT_COLUMNS)) + "\n"
                )
                continue

            try:
                statistics = stats_over_range(cpg_ranges, chrom, start, end)
            except Exception as exc:
                raise BetaStatsError(
                    f"Failed to summarize region on line {line_number}: {exc}"
                ) from exc

            formatted = [format_stat(value, na_rep) for value in statistics]
            if len(formatted) != len(STAT_COLUMNS):
                raise BetaStatsError(
                    "stats_over_range() returned "
                    f"{len(formatted)} values; expected {len(STAT_COLUMNS)}."
                )

            output_handle.write(line + "\t" + "\t".join(formatted) + "\n")
            processed += 1

    return processed, invalid


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_file = Path(f"{args.out_prefix}.region_stats.tsv")

    try:
        printlog(f'Reading CpG methylation file: "{args.input_file}"')
        cpg_ranges = read_CpG_bed(str(args.input_file))

        printlog(f'Reading genomic regions: "{args.region}"')
        printlog(f'Writing regional methylation statistics: "{output_file}"')

        processed, invalid = summarize_regions(
            cpg_ranges=cpg_ranges,
            region_file=args.region,
            output_file=output_file,
            write_header=args.header,
            na_rep=args.na_rep,
        )

        printlog(f"Processed {processed} valid genomic regions.")
        if invalid:
            printlog(
                f"Warning: {invalid} invalid region lines were written with "
                f"{args.na_rep!r} statistics."
            )

    except BetaStatsError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

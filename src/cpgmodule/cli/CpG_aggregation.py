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

"""Aggregate CpG methylation values over user-defined genomic regions.

Two input data types are supported in column 4 of the CpG BED file:

* ``count``: methylated,total read counts, e.g. ``3,10``
* ``beta``: methylation beta values, e.g. ``0.30``

For count data, CpG-level outliers can optionally be removed using a
Bonferroni-adjusted two-sided exact binomial test. The null methylation
probability for each region is estimated from all CpGs overlapping that region.

The region file may contain an arbitrary number of columns. All original region
columns are preserved in the output, and aggregation columns are appended.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
from bx.intervals import Intersecter, Interval
from scipy.stats import binomtest

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class AggregationError(RuntimeError):
    """Raised when input validation or CpG aggregation fails."""


COUNT_COLUMNS = (
    "N_CpG_filtered",
    "N_methyl_filtered",
    "N_total_filtered",
    "Beta_filtered",
    "N_CpG_original",
    "N_methyl_original",
    "N_total_original",
    "Beta_original",
)

BETA_COLUMNS = (
    "N_CpG",
    "Mean_beta",
    "Median_beta",
    "Min_beta",
    "Max_beta",
    "Std_beta",
)


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_aggregation",
        description=(
            "Aggregate CpG methylation values over user-defined genomic regions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "CpG BED3+ file. Column 4 contains either count values such as "
            "'3,10' or beta values such as '0.30'. Compressed input is supported."
        ),
    )
    parser.add_argument(
        "-b",
        "--bed",
        dest="region_file",
        required=True,
        type=Path,
        help=(
            "BED3+ file defining genomic regions. The file may contain an "
            "arbitrary number of columns; all original columns are preserved."
        ),
    )
    parser.add_argument(
        "-t",
        "--type",
        dest="data_type",
        required=True,
        choices=("count", "beta"),
        help="Data type stored in column 4 of the CpG input file.",
    )
    parser.add_argument(
        "-a",
        "--alpha",
        dest="alpha",
        type=float,
        default=0.05,
        help=(
            "Family-wise error rate used for CpG outlier filtering in count "
            "mode. The per-CpG cutoff is alpha / number_of_CpGs_in_region."
        ),
    )
    parser.add_argument(
        "--no_outlier_filter",
        action="store_true",
        help="Disable CpG outlier filtering in count mode.",
    )
    parser.add_argument(
        "--min_cpg",
        type=int,
        default=1,
        help="Minimum number of overlapping valid CpGs required for a region.",
    )
    parser.add_argument(
        "--ddof",
        type=int,
        choices=(0, 1),
        default=0,
        help="Delta degrees of freedom for beta-value standard deviation.",
    )
    parser.add_argument(
        "--header",
        action="store_true",
        help=(
            "Write a header. Original region columns are named chrom/start/end "
            "followed by column_4, column_5, ... when no explicit names exist."
        ),
    )
    parser.add_argument(
        "--na_rep",
        default="NA",
        help="Text used for unavailable output values.",
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
        parser.error(f"Input CpG file does not exist: {args.input_file}")
    if not args.region_file.is_file():
        parser.error(f"Region BED file does not exist: {args.region_file}")
    if not 0.0 <= args.alpha <= 1.0:
        parser.error("--alpha must be between 0 and 1")
    if args.min_cpg < 1:
        parser.error("--min_cpg must be at least 1")
    if "\t" in args.na_rep or "\n" in args.na_rep:
        parser.error("--na_rep cannot contain tabs or newlines")


def iter_data_lines(path: Path):
    """Yield non-empty, non-directive lines from a BED-like file."""
    for line_number, raw_line in enumerate(ireader.reader(str(path)), start=1):
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(("#", "track", "browser")):
            continue
        yield line_number, line


def parse_count(value: str, line_number: int) -> tuple[int, int]:
    """Parse a methylated,total count pair."""
    try:
        methyl_text, total_text = value.split(",", 1)
        methylated = int(methyl_text)
        total = int(total_text)
    except Exception as exc:
        raise AggregationError(
            f"Invalid count value {value!r} on line {line_number}; "
            "expected 'methylated,total', e.g. '3,10'."
        ) from exc

    if methylated < 0 or total < 0:
        raise AggregationError(
            f"Negative count found on line {line_number}: {value!r}"
        )
    if methylated > total:
        raise AggregationError(
            f"Methylated reads exceed total reads on line {line_number}: {value!r}"
        )
    return methylated, total


def parse_beta(value: str, line_number: int) -> float:
    """Parse and validate a beta value."""
    try:
        beta = float(value)
    except ValueError as exc:
        raise AggregationError(
            f"Invalid beta value {value!r} on line {line_number}."
        ) from exc

    if not math.isfinite(beta) or not 0.0 <= beta <= 1.0:
        raise AggregationError(
            f"Beta value outside [0, 1] on line {line_number}: {value!r}"
        )
    return beta


def build_interval_tree(path: Path, data_type: str) -> dict[str, Intersecter]:
    """Build chromosome-specific interval trees from the CpG BED file."""
    printlog(f'Reading CpG file: "{path}"')
    trees: dict[str, Intersecter] = {}
    valid_records = 0

    for line_number, line in iter_data_lines(path):
        fields = line.split()
        if len(fields) < 4:
            raise AggregationError(
                f"Line {line_number} of {path} has fewer than four columns."
            )

        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise AggregationError(
                f"Invalid BED coordinates on line {line_number}: {line}"
            ) from exc

        if start < 0 or end <= start:
            raise AggregationError(
                f"Invalid BED interval on line {line_number}: {chrom}:{start}-{end}"
            )

        if data_type == "count":
            value = parse_count(fields[3], line_number)
        else:
            value = parse_beta(fields[3], line_number)

        tree = trees.setdefault(chrom, Intersecter())
        tree.add_interval(Interval(start, end, value=value))
        valid_records += 1

    if valid_records == 0:
        raise AggregationError("No valid CpG records were found in the input file.")

    printlog(f"Loaded {valid_records} CpG records.")
    return trees


def find_overlaps(
    chrom: str,
    start: int,
    end: int,
    trees: dict[str, Intersecter],
) -> list:
    """Return CpG values overlapping one genomic region."""
    if chrom not in trees:
        return []
    return [interval.value for interval in trees[chrom].find(start, end)]


def aggregate_counts(
    hits: list[tuple[int, int]],
    alpha: float,
    filter_outliers: bool,
    min_cpg: int,
) -> Optional[list[float | int]]:
    """Aggregate count data and optionally remove CpG-level outliers."""
    valid_hits = [(m, t) for m, t in hits if t > 0]
    if len(valid_hits) < min_cpg:
        return None

    methyl = np.asarray([m for m, _ in valid_hits], dtype=int)
    total = np.asarray([t for _, t in valid_hits], dtype=int)

    original_cpg = int(len(total))
    original_methyl = int(methyl.sum())
    original_total = int(total.sum())

    if original_total == 0:
        return None

    original_beta = original_methyl / original_total

    keep = np.ones(original_cpg, dtype=bool)

    if filter_outliers and original_cpg > 1 and 0.0 < original_beta < 1.0:
        p_cutoff = alpha / original_cpg

        for index, (methylated, coverage) in enumerate(valid_hits):
            # Exact two-sided binomial test. This is preferable to comparing
            # binom.cdf() and 1-cdf(), which does not correctly represent the
            # upper tail for discrete data.
            p_value = binomtest(
                k=methylated,
                n=coverage,
                p=original_beta,
                alternative="two-sided",
            ).pvalue
            if p_value < p_cutoff:
                keep[index] = False

    filtered_methyl = methyl[keep]
    filtered_total = total[keep]

    filtered_cpg = int(len(filtered_total))
    if filtered_cpg < min_cpg or filtered_total.sum() == 0:
        filtered_methyl_sum = 0
        filtered_total_sum = 0
        filtered_beta = np.nan
    else:
        filtered_methyl_sum = int(filtered_methyl.sum())
        filtered_total_sum = int(filtered_total.sum())
        filtered_beta = filtered_methyl_sum / filtered_total_sum

    return [
        filtered_cpg,
        filtered_methyl_sum,
        filtered_total_sum,
        filtered_beta,
        original_cpg,
        original_methyl,
        original_total,
        original_beta,
    ]


def aggregate_betas(
    hits: list[float],
    min_cpg: int,
    ddof: int,
) -> Optional[list[float | int]]:
    """Aggregate beta values over one region."""
    values = np.asarray(hits, dtype=float)
    values = values[np.isfinite(values)]

    if values.size < min_cpg:
        return None

    std_value = np.std(values, ddof=ddof) if values.size > ddof else np.nan

    return [
        int(values.size),
        float(np.mean(values)),
        float(np.median(values)),
        float(np.min(values)),
        float(np.max(values)),
        float(std_value),
    ]


def infer_region_column_count(region_file: Path) -> int:
    """Infer the number of columns from the first valid region record."""
    for _, line in iter_data_lines(region_file):
        return len(line.split())
    raise AggregationError("Region file contains no valid data records.")


def build_region_header(region_file: Path, stat_columns: Sequence[str]) -> list[str]:
    """Create a generic header for arbitrary-column BED-like input."""
    n_columns = infer_region_column_count(region_file)
    if n_columns < 3:
        raise AggregationError("Region file must contain at least three columns.")

    columns = ["chrom", "start", "end"]
    columns.extend(f"column_{i}" for i in range(4, n_columns + 1))
    columns.extend(stat_columns)
    return columns


def format_value(value, na_rep: str) -> str:
    """Format one output value."""
    if value is None:
        return na_rep
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)

    if not math.isfinite(numeric):
        return na_rep
    if numeric.is_integer():
        return str(int(numeric))
    return f"{numeric:.6f}"


def process_regions(
    region_file: Path,
    trees: dict[str, Intersecter],
    data_type: str,
    output_file: Path,
    *,
    alpha: float,
    filter_outliers: bool,
    min_cpg: int,
    ddof: int,
    write_header: bool,
    na_rep: str,
) -> tuple[int, int]:
    """Aggregate CpGs for all target regions and write output."""
    stat_columns = COUNT_COLUMNS if data_type == "count" else BETA_COLUMNS
    expected_columns = infer_region_column_count(region_file)
    processed = 0
    invalid = 0

    with output_file.open("w", encoding="utf-8") as handle:
        if write_header:
            handle.write(
                "\t".join(build_region_header(region_file, stat_columns)) + "\n"
            )

        for line_number, line in iter_data_lines(region_file):
            fields = line.split()

            if len(fields) != expected_columns:
                invalid += 1
                handle.write(
                    line
                    + "\t"
                    + "\t".join([na_rep] * len(stat_columns))
                    + "\n"
                )
                continue

            if len(fields) < 3:
                invalid += 1
                handle.write(
                    line
                    + "\t"
                    + "\t".join([na_rep] * len(stat_columns))
                    + "\n"
                )
                continue

            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                invalid += 1
                handle.write(
                    line
                    + "\t"
                    + "\t".join([na_rep] * len(stat_columns))
                    + "\n"
                )
                continue

            if start < 0 or end <= start:
                invalid += 1
                handle.write(
                    line
                    + "\t"
                    + "\t".join([na_rep] * len(stat_columns))
                    + "\n"
                )
                continue

            hits = find_overlaps(chrom, start, end, trees)

            if data_type == "count":
                result = aggregate_counts(
                    hits,
                    alpha=alpha,
                    filter_outliers=filter_outliers,
                    min_cpg=min_cpg,
                )
            else:
                result = aggregate_betas(
                    hits,
                    min_cpg=min_cpg,
                    ddof=ddof,
                )

            if result is None:
                output_values = [na_rep] * len(stat_columns)
            else:
                output_values = [format_value(value, na_rep) for value in result]

            handle.write(line + "\t" + "\t".join(output_values) + "\n")
            processed += 1

    return processed, invalid


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_file = Path(f"{args.out_prefix}.aggregation.tsv")

    try:
        trees = build_interval_tree(args.input_file, args.data_type)

        printlog(f'Reading target regions: "{args.region_file}"')
        printlog(f'Writing aggregated methylation data: "{output_file}"')

        processed, invalid = process_regions(
            region_file=args.region_file,
            trees=trees,
            data_type=args.data_type,
            output_file=output_file,
            alpha=args.alpha,
            filter_outliers=not args.no_outlier_filter,
            min_cpg=args.min_cpg,
            ddof=args.ddof,
            write_header=args.header,
            na_rep=args.na_rep,
        )

        printlog(f"Processed {processed} genomic regions.")
        if invalid:
            printlog(
                f"Warning: {invalid} region lines were invalid or had an "
                "inconsistent number of columns and were written with missing "
                "summary values."
            )

    except AggregationError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

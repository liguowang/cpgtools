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

"""Annotate CpG sites by overlap with user-defined genomic regions.

The CpG input file and annotation file must each contain at least three BED
columns: chromosome, start, and end. Additional columns are preserved.

If the annotation file contains a fourth column, that value is used as the
annotation label. Otherwise, labels are generated from genomic coordinates
(e.g. ``chr1:100-200``).

When multiple annotation regions overlap the same CpG, unique labels are
combined using a configurable separator.

A centered window can optionally be used for each annotation region. For
example, ``--window 1000`` restricts each annotation interval to a 1000-bp
window centered on the original region midpoint, clipped to the original
annotation interval.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from bx.intervals import Intersecter, Interval

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class AnnotationError(RuntimeError):
    """Raised when annotation input or interval processing fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_anno_position",
        description="Annotate CpG sites by overlap with BED-defined genomic regions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input CpG BED3+ file. Additional columns are preserved. "
            "Compressed input is supported by the CpGtools reader."
        ),
    )
    parser.add_argument(
        "-a",
        "--annotation",
        dest="annotation_file",
        required=True,
        type=Path,
        help=(
            "Annotation BED3+ file. Column 4 is used as the annotation label "
            "when present; otherwise chromosome:start-end is used."
        ),
    )
    parser.add_argument(
        "-w",
        "--window",
        dest="window_size",
        type=int,
        default=0,
        help=(
            "Size in bp of a window centered on each annotation-region midpoint. "
            "A value of 0 uses the full annotation interval."
        ),
    )
    parser.add_argument(
        "--separator",
        default=",",
        help="Separator used to join multiple overlapping annotation labels.",
    )
    parser.add_argument(
        "--na_rep",
        default="NA",
        help="Text written when a CpG has no overlapping annotation.",
    )
    parser.add_argument(
        "--annotation_name",
        default=None,
        help=(
            "Name of the appended annotation column when --header is used. "
            "Default: annotation filename."
        ),
    )
    parser.add_argument(
        "--header",
        "-l",
        action="store_true",
        help="Treat the first non-directive line of the CpG input as a header.",
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

    if not args.annotation_file.is_file():
        parser.error(f"Annotation BED file does not exist: {args.annotation_file}")

    if args.window_size < 0:
        parser.error("--window must be non-negative")

    if not args.separator:
        parser.error("--separator cannot be empty")

    if "\t" in args.separator or "\n" in args.separator:
        parser.error("--separator cannot contain tabs or newlines")

    if "\t" in args.na_rep or "\n" in args.na_rep:
        parser.error("--na_rep cannot contain tabs or newlines")


def iter_data_lines(path: Path):
    """Yield BED-like data lines while skipping comments and UCSC directives."""
    for line_number, raw_line in enumerate(ireader.reader(str(path)), start=1):
        line = raw_line.rstrip("\r\n")
        stripped = line.strip()

        if not stripped:
            continue
        if stripped.startswith(("#", "track", "browser")):
            continue

        yield line_number, line


def centered_interval(start: int, end: int, window_size: int) -> tuple[int, int]:
    """Return a centered interval clipped to the original BED interval."""
    if window_size <= 0 or window_size >= (end - start):
        return start, end

    midpoint = start + (end - start) // 2

    left = window_size // 2
    right = window_size - left

    window_start = max(start, midpoint - left)
    window_end = min(end, midpoint + right)

    if window_end <= window_start:
        return start, end

    return window_start, window_end


def build_interval_tree(
    annotation_file: Path,
    window_size: int = 0,
) -> tuple[dict[str, Intersecter], int]:
    """Build chromosome-specific interval trees from an annotation BED file."""
    printlog(f'Building interval tree from annotation file: "{annotation_file}"')

    trees: dict[str, Intersecter] = {}
    loaded = 0

    for line_number, line in iter_data_lines(annotation_file):
        fields = line.split()

        if len(fields) < 3:
            raise AnnotationError(
                f"Annotation line {line_number} has fewer than three columns."
            )

        chrom = fields[0]

        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise AnnotationError(
                f"Invalid annotation coordinates on line {line_number}: {line}"
            ) from exc

        if start < 0 or end <= start:
            raise AnnotationError(
                f"Invalid annotation interval on line {line_number}: "
                f"{chrom}:{start}-{end}"
            )

        label = fields[3] if len(fields) >= 4 else f"{chrom}:{start}-{end}"

        interval_start, interval_end = centered_interval(
            start,
            end,
            window_size,
        )

        tree = trees.setdefault(chrom, Intersecter())
        tree.add_interval(
            Interval(
                interval_start,
                interval_end,
                value=label,
            )
        )
        loaded += 1

    if loaded == 0:
        raise AnnotationError(
            f"No valid annotation regions were found in {annotation_file}."
        )

    printlog(f"Loaded {loaded} annotation regions.")
    return trees, loaded


def find_intervals(
    chrom: str,
    start: int,
    end: int,
    trees: dict[str, Intersecter],
) -> list[str]:
    """Return sorted unique labels overlapping one CpG interval."""
    if chrom not in trees:
        return []

    overlaps = trees[chrom].find(start, end)
    return sorted({str(interval.value) for interval in overlaps})


def annotate_cpg_file(
    input_file: Path,
    trees: dict[str, Intersecter],
    output_file: Path,
    *,
    has_header: bool,
    annotation_column_name: str,
    separator: str,
    na_rep: str,
) -> tuple[int, int, int]:
    """Annotate all CpG records and write the resulting table."""
    processed = 0
    annotated = 0
    invalid = 0
    first_data_line = True

    with output_file.open("w", encoding="utf-8") as output_handle:
        for line_number, line in iter_data_lines(input_file):
            fields = line.split()

            if first_data_line and has_header:
                output_handle.write(
                    line + "\t" + annotation_column_name + "\n"
                )
                first_data_line = False
                continue

            first_data_line = False

            if len(fields) < 3:
                invalid += 1
                output_handle.write(line + "\t" + na_rep + "\n")
                continue

            chrom = fields[0]

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                invalid += 1
                output_handle.write(line + "\t" + na_rep + "\n")
                continue

            if start < 0 or end <= start:
                invalid += 1
                output_handle.write(line + "\t" + na_rep + "\n")
                continue

            overlaps = find_intervals(
                chrom=chrom,
                start=start,
                end=end,
                trees=trees,
            )

            if overlaps:
                annotation = separator.join(overlaps)
                annotated += 1
            else:
                annotation = na_rep

            output_handle.write(line + "\t" + annotation + "\n")
            processed += 1

    return processed, annotated, invalid


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_file = Path(f"{args.out_prefix}.anno.tsv")

    annotation_column_name = (
        args.annotation_name
        if args.annotation_name
        else args.annotation_file.name
    )

    try:
        trees, _ = build_interval_tree(
            annotation_file=args.annotation_file,
            window_size=args.window_size,
        )

        printlog(f'Reading CpG file: "{args.input_file}"')
        printlog(f'Writing annotated CpGs: "{output_file}"')

        processed, annotated, invalid = annotate_cpg_file(
            input_file=args.input_file,
            trees=trees,
            output_file=output_file,
            has_header=args.header,
            annotation_column_name=annotation_column_name,
            separator=args.separator,
            na_rep=args.na_rep,
        )

        printlog(f"Processed CpGs: {processed}")
        printlog(f"CpGs with annotation: {annotated}")
        printlog(f"CpGs without annotation: {processed - annotated}")

        if invalid:
            printlog(
                f"Warning: {invalid} invalid CpG lines were written with "
                f"{args.na_rep!r} annotation."
            )

    except AnnotationError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

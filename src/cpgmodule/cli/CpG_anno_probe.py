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

"""Append array-probe annotations to a tabular input file.

The input file may contain arbitrary columns, with one column holding Illumina
DNA methylation probe IDs (for example, cg00000029). The annotation file must
contain a probe-ID column followed by one or more annotation columns.

The original input columns are preserved and annotation columns are appended.
Compressed input files are supported through the CpGtools reader.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class ProbeAnnotationError(RuntimeError):
    """Raised when probe annotation input or processing fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_anno_probe",
        description="Annotate CpGs by their probe ID (using a prebuilt tabular file).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input tabular file containing probe IDs in one column. "
            "Compressed input is supported."
        ),
    )
    parser.add_argument(
        "-a",
        "--annotation",
        dest="annotation_file",
        required=True,
        type=Path,
        help=(
            "Probe annotation file. The first column must contain probe IDs; "
            "remaining columns contain annotation fields."
        ),
    )
    parser.add_argument(
        "-p",
        "--probe_column",
        dest="probe_column",
        type=int,
        default=0,
        help="Zero-based column index containing probe IDs in the input file.",
    )
    parser.add_argument(
        "--annotation_probe_column",
        type=int,
        default=0,
        help="Zero-based probe-ID column index in the annotation file.",
    )
    parser.add_argument(
        "--header",
        "-l",
        action="store_true",
        help="Treat the first line of the input file as a header.",
    )
    parser.add_argument(
        "--annotation_header",
        choices=("auto", "yes", "no"),
        default="auto",
        help=(
            "Whether the annotation file contains a header. 'auto' detects "
            "common probe-ID header names such as probeID or IlmnID."
        ),
    )
    parser.add_argument(
        "--na_rep",
        default="NA",
        help="Text written when a probe has no matching annotation.",
    )
    parser.add_argument(
        "--duplicate_policy",
        choices=("first", "last", "error"),
        default="error",
        help="How to handle duplicate probe IDs in the annotation file.",
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

    if not args.annotation_file.is_file():
        parser.error(f"Annotation file does not exist: {args.annotation_file}")

    if args.probe_column < 0:
        parser.error("--probe_column must be non-negative")

    if args.annotation_probe_column < 0:
        parser.error("--annotation_probe_column must be non-negative")

    if "\t" in args.na_rep or "\n" in args.na_rep:
        parser.error("--na_rep cannot contain tabs or newlines")


def split_fields(line: str) -> list[str]:
    """Split a tabular line while preferring tabs over generic whitespace."""
    stripped = line.rstrip("\r\n")
    if "\t" in stripped:
        return stripped.split("\t")
    return stripped.split()


def detect_annotation_header(fields: list[str], probe_col: int) -> bool:
    """Detect common annotation probe-ID headers."""
    if probe_col >= len(fields):
        return False

    token = fields[probe_col].strip().lower()
    return token in {
        "probeid",
        "probe_id",
        "ilmnid",
        "ilmn_id",
        "cpg",
        "cpg_id",
        "cpgid",
        "name",
    }


def read_annotation(
    annotation_file: Path,
    probe_column: int,
    header_mode: str,
    duplicate_policy: str,
) -> tuple[list[str], dict[str, list[str]], int]:
    """Read probe annotations and return headers and probe-indexed values."""
    printlog(f'Reading annotation file: "{annotation_file}"')

    annotation_data: dict[str, list[str]] = {}
    annotation_headers: list[str] = []
    duplicate_count = 0
    first_record = True
    expected_columns: Optional[int] = None

    for line_number, raw_line in enumerate(
        ireader.reader(str(annotation_file)),
        start=1,
    ):
        line = raw_line.rstrip("\r\n")

        if not line.strip():
            continue
        if line.lstrip().startswith(("#", "track", "browser")):
            continue

        fields = split_fields(line)

        if expected_columns is None:
            expected_columns = len(fields)

        if len(fields) != expected_columns:
            raise ProbeAnnotationError(
                f"Annotation file has inconsistent column counts: expected "
                f"{expected_columns}, found {len(fields)} on line {line_number}."
            )

        if probe_column >= len(fields):
            raise ProbeAnnotationError(
                f"--annotation_probe_column={probe_column} exceeds the "
                f"{len(fields)} columns on annotation line {line_number}."
            )

        if first_record:
            has_header = (
                header_mode == "yes"
                or (
                    header_mode == "auto"
                    and detect_annotation_header(fields, probe_column)
                )
            )

            if has_header:
                annotation_headers = [
                    value
                    for index, value in enumerate(fields)
                    if index != probe_column
                ]
                first_record = False
                continue

            annotation_headers = [
                f"annotation_{index}"
                for index in range(1, len(fields))
            ]
            first_record = False

        probe_id = fields[probe_column].strip()
        if not probe_id:
            continue

        values = [
            value
            for index, value in enumerate(fields)
            if index != probe_column
        ]

        if probe_id in annotation_data:
            duplicate_count += 1

            if duplicate_policy == "error":
                raise ProbeAnnotationError(
                    f"Duplicate probe ID {probe_id!r} found in annotation "
                    f"file on line {line_number}."
                )

            if duplicate_policy == "first":
                continue

        annotation_data[probe_id] = values

    if expected_columns is None:
        raise ProbeAnnotationError("Annotation file contains no data records.")

    if not annotation_headers:
        annotation_headers = [
            f"annotation_{index}"
            for index in range(1, expected_columns)
        ]

    printlog(f"Loaded annotations for {len(annotation_data)} probes.")
    if duplicate_count:
        printlog(
            f"Encountered {duplicate_count} duplicate annotation probe IDs "
            f"using policy {duplicate_policy!r}."
        )

    return annotation_headers, annotation_data, duplicate_count


def annotate_input(
    input_file: Path,
    output_file: Path,
    *,
    probe_column: int,
    has_header: bool,
    annotation_headers: list[str],
    annotation_data: dict[str, list[str]],
    na_rep: str,
) -> tuple[int, int, int]:
    """Append annotations to each input record."""
    processed = 0
    matched = 0
    invalid = 0
    first_record = True
    expected_columns: Optional[int] = None

    missing_annotation = [na_rep] * len(annotation_headers)

    with output_file.open("w", encoding="utf-8") as output_handle:
        for line_number, raw_line in enumerate(
            ireader.reader(str(input_file)),
            start=1,
        ):
            line = raw_line.rstrip("\r\n")

            if not line.strip():
                continue

            fields = split_fields(line)

            if first_record and has_header:
                expected_columns = len(fields)
                output_handle.write(
                    "\t".join(fields + annotation_headers) + "\n"
                )
                first_record = False
                continue

            first_record = False

            if expected_columns is None:
                expected_columns = len(fields)

            if len(fields) != expected_columns:
                invalid += 1
                output_handle.write(
                    "\t".join(fields + missing_annotation) + "\n"
                )
                continue

            if probe_column >= len(fields):
                raise ProbeAnnotationError(
                    f"--probe_column={probe_column} exceeds the {len(fields)} "
                    f"columns on input line {line_number}."
                )

            probe_id = fields[probe_column].strip()
            annotation = annotation_data.get(probe_id)

            if annotation is None:
                annotation = missing_annotation
            else:
                matched += 1

            output_handle.write(
                "\t".join(fields + annotation) + "\n"
            )
            processed += 1

    return processed, matched, invalid


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_file = Path(f"{args.out_prefix}.anno.tsv")

    try:
        (
            annotation_headers,
            annotation_data,
            _,
        ) = read_annotation(
            annotation_file=args.annotation_file,
            probe_column=args.annotation_probe_column,
            header_mode=args.annotation_header,
            duplicate_policy=args.duplicate_policy,
        )

        printlog(f'Annotating input file: "{args.input_file}"')
        printlog(f'Writing annotated table: "{output_file}"')

        processed, matched, invalid = annotate_input(
            input_file=args.input_file,
            output_file=output_file,
            probe_column=args.probe_column,
            has_header=args.header,
            annotation_headers=annotation_headers,
            annotation_data=annotation_data,
            na_rep=args.na_rep,
        )

        printlog(f"Processed records: {processed}")
        printlog(f"Matched probes: {matched}")
        printlog(f"Unmatched probes: {processed - matched}")

        if invalid:
            printlog(
                f"Warning: {invalid} input lines had inconsistent column "
                f"counts and were written with {args.na_rep!r} annotations."
            )

    except ProbeAnnotationError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

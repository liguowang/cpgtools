#!/usr/bin/env python3
#
# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: wangliguo78@gmail.com
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

"""
Description
-----------
Perform differential CpG analysis using Fisher's exact test.

This command is intended for two-group comparison of methylated and total read
counts. If a group contains multiple replicates, methylated and total reads are
summed across replicates before Fisher's exact test is applied.

Input count values must be written as:

    methylated_count,total_count

where both values are non-negative integers and
``methylated_count <= total_count``.

Example
-------
    cgID    sample_1   sample_2
    CpG_1   129,170    166,178
    CpG_2   24,77      67,99

Output
------
Three columns are appended to the original input table:

* OddsRatio
* pval
* adj.pval

Multiple-testing correction uses the existing CpGtools
``padjust.multiple_testing_correction`` implementation.
"""

import argparse
import collections
import re
import sys
from pathlib import Path

from scipy import stats

from cpgmodule import ireader
from cpgmodule import padjust
from cpgmodule._version import __version__
from cpgmodule.utils import printlog, read_grp_file1


COUNT_RE = re.compile(r"^(\d+)\s*,\s*(\d+)$")


def build_parser():
    """Build and return the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help=(
            "Data file containing methylation counts as "
            "'methylated_count,total_count'. The first row contains sample "
            "IDs and the first column contains CpG/probe IDs. Plain text, "
            ".gz, .bz2, and other formats supported by cpgmodule.ireader "
            "may be used."
        ),
    )

    parser.add_argument(
        "-g",
        "--group",
        dest="group_file",
        required=True,
        help=(
            "Comma-separated two-column group file with a header. "
            "Column 1 contains sample IDs and column 2 contains group IDs. "
            "Exactly two groups are required."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help="Output prefix. Results are written to <prefix>.pval.txt.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def parse_count(value):
    """
    Parse one methylated_count,total_count value.

    Returns
    -------
    tuple[int, int] or None
        Parsed counts, or None if the field is not a valid count expression.

    Raises
    ------
    ValueError
        If the value is syntactically valid but biologically impossible.
    """
    match = COUNT_RE.match(value)

    if match is None:
        return None

    methylated = int(match.group(1))
    total = int(match.group(2))

    if total <= 0:
        raise ValueError(
            f"total_count must be greater than zero: {value}"
        )

    if methylated > total:
        raise ValueError(
            f"methylated_count cannot exceed total_count: {value}"
        )

    return methylated, total


def read_groups(group_file):
    """
    Read and validate the two-group sample mapping.

    Returns
    -------
    tuple
        sample-to-group mapping, ordered group IDs, group-to-samples mapping.
    """
    printlog(f'Read group file "{group_file}" ...')

    samples, groups = read_grp_file1(group_file)

    if not samples:
        raise ValueError("Group file contains no samples")

    sample_to_group = dict(zip(samples, groups))

    if len(sample_to_group) != len(samples):
        raise ValueError("Sample IDs in group file are not unique")

    group_to_samples = collections.defaultdict(list)
    for sample, group in zip(samples, groups):
        group_to_samples[group].append(sample)

    group_ids = sorted(group_to_samples)

    for group_id in group_ids:
        print(
            f"\tGroup {group_id} has "
            f"{len(group_to_samples[group_id])} samples:",
            file=sys.stderr,
        )
        print(
            "\t\t" + ",".join(group_to_samples[group_id]),
            file=sys.stderr,
        )

    if len(group_ids) != 2:
        raise ValueError("Exactly two groups are required")

    return sample_to_group, group_ids, group_to_samples


def analyze_probe(
    proportions,
    sample_ids,
    sample_to_group,
    group_ids,
    line_number,
    probe_id,
):
    """
    Aggregate replicate counts by group and run Fisher's exact test.

    Returns
    -------
    tuple[float, float]
        Odds ratio and p-value.
    """
    group_counts = {
        group_id: {"methyl": 0, "unmethyl": 0}
        for group_id in group_ids
    }

    for sample_id, value in zip(sample_ids, proportions):
        group_id = sample_to_group.get(sample_id)

        # Samples not defined in the group file are ignored.
        if group_id is None:
            continue

        try:
            parsed = parse_count(value)
        except ValueError as exc:
            raise ValueError(
                f'Invalid count data on line {line_number}, '
                f'probe "{probe_id}", sample "{sample_id}": {exc}'
            ) from exc

        # Preserve the original behavior: non-count fields are treated as
        # missing and ignored.
        if parsed is None:
            continue

        methylated, total = parsed

        group_counts[group_id]["methyl"] += methylated
        group_counts[group_id]["unmethyl"] += total - methylated

    table = [
        [
            group_counts[group_ids[0]]["methyl"],
            group_counts[group_ids[0]]["unmethyl"],
        ],
        [
            group_counts[group_ids[1]]["methyl"],
            group_counts[group_ids[1]]["unmethyl"],
        ],
    ]

    # scipy.stats.fisher_exact returns NaN when the contingency table is
    # completely empty in newer SciPy versions. Explicitly preserve this as
    # a missing statistical result rather than failing.
    odds_ratio, p_value = stats.fisher_exact(table)

    return float(odds_ratio), float(p_value)


def analyze_file(input_file, sample_to_group, group_ids):
    """
    Analyze all CpG/probe rows in the input matrix.

    Returns
    -------
    tuple
        Original header, original data lines, probe IDs, odds ratios, p-values.
    """
    header = None
    sample_ids = None
    records = []
    probe_ids = []
    odds_ratios = []
    p_values = []

    for line_number, raw_line in enumerate(
        ireader.reader(input_file),
        start=1,
    ):
        line = raw_line.rstrip("\r\n")

        if not line:
            continue

        fields = line.split()

        if header is None:
            if len(fields) < 2:
                raise ValueError(
                    "Input matrix header must contain at least one sample ID"
                )

            header = line
            sample_ids = fields[1:]

            if len(sample_ids) != len(set(sample_ids)):
                raise ValueError(
                    "Sample IDs in input matrix are not unique"
                )

            missing = [
                sample
                for sample in sample_to_group
                if sample not in sample_ids
            ]
            if missing:
                raise ValueError(
                    "Cannot find sample ID(s) in data file: "
                    + ", ".join(missing)
                )

            continue

        probe_id = fields[0]
        proportions = fields[1:]

        if len(proportions) != len(sample_ids):
            print(
                f"Warning: line {line_number} has {len(proportions)} data "
                f"fields but header has {len(sample_ids)} samples. "
                "Missing fields are treated as missing; extra fields are "
                "ignored.",
                file=sys.stderr,
            )

        if len(proportions) < len(sample_ids):
            proportions = proportions + ["NA"] * (
                len(sample_ids) - len(proportions)
            )
        elif len(proportions) > len(sample_ids):
            proportions = proportions[: len(sample_ids)]

        odds_ratio, p_value = analyze_probe(
            proportions=proportions,
            sample_ids=sample_ids,
            sample_to_group=sample_to_group,
            group_ids=group_ids,
            line_number=line_number,
            probe_id=probe_id,
        )

        records.append(line)
        probe_ids.append(probe_id)
        odds_ratios.append(odds_ratio)
        p_values.append(p_value)

    if header is None:
        raise ValueError("Input data file is empty")

    if not records:
        raise ValueError("Input data file contains no CpG/probe rows")

    return header, records, probe_ids, odds_ratios, p_values


def adjust_pvalues(p_values):
    """
    Apply CpGtools multiple-testing correction while preserving NaNs.

    Returns
    -------
    list[float]
        Adjusted p-values in original order.
    """
    valid_indices = []
    valid_pvalues = []

    for index, p_value in enumerate(p_values):
        if p_value == p_value:  # NaN-safe finite/missing check
            valid_indices.append(index)
            valid_pvalues.append(p_value)

    adjusted = [float("nan")] * len(p_values)

    if valid_pvalues:
        corrected = padjust.multiple_testing_correction(valid_pvalues)
        for index, q_value in zip(valid_indices, corrected):
            adjusted[index] = float(q_value)

    return adjusted


def format_stat(value):
    """Format a numeric statistic while preserving missing values."""
    if value != value:
        return "NA"
    return str(value)


def run_fisher(input_file, group_file, output_prefix):
    """
    Run Fisher's exact test differential methylation analysis.

    Returns
    -------
    str
        Path to the generated output file.
    """
    sample_to_group, group_ids, _ = read_groups(group_file)

    (
        header,
        records,
        probe_ids,
        odds_ratios,
        p_values,
    ) = analyze_file(
        input_file=input_file,
        sample_to_group=sample_to_group,
        group_ids=group_ids,
    )

    printlog("Perform Benjamini-Hochberg (FDR) correction ...")
    adjusted_pvalues = adjust_pvalues(p_values)

    output_path = Path(f"{output_prefix}.pval.txt")
    if output_path.parent != Path("."):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    printlog(f'Writing to "{output_path}"')

    with open(output_path, "w") as fout:
        print(
            header + "\tOddsRatio\tpval\tadj.pval",
            file=fout,
        )

        for (
            line,
            probe_id,
            odds_ratio,
            p_value,
            adjusted_pvalue,
        ) in zip(
            records,
            probe_ids,
            odds_ratios,
            p_values,
            adjusted_pvalues,
        ):
            print(
                line
                + "\t"
                + format_stat(odds_ratio)
                + "\t"
                + format_stat(p_value)
                + "\t"
                + format_stat(adjusted_pvalue),
                file=fout,
            )

    return str(output_path)


def main(argv=None):
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str] or None
        Optional argument list. When None, argparse reads ``sys.argv``.
        Passing a list allows this function to be called directly from the
        CpGtools dispatcher or from tests.

    Returns
    -------
    int
        Process-style return code.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        run_fisher(
            input_file=args.input_file,
            group_file=args.group_file,
            output_prefix=args.out_file,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

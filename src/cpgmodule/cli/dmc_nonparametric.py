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
Perform differential CpG analysis using nonparametric tests.

For exactly two groups, the Mann-Whitney U test is used.
For three or more groups, the Kruskal-Wallis H-test is used.

The first row of the input matrix must contain unique sample IDs, and the
first column must contain CpG/probe IDs. Non-numeric values in the data matrix
are treated as missing and ignored.

Output
------
Two columns are appended to the original input table:

* pval
* adj.pval

Multiple-testing correction uses the existing CpGtools
``padjust.multiple_testing_correction`` implementation.
"""

import argparse
import collections
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

from cpgmodule import ireader
from cpgmodule import padjust
from cpgmodule._version import __version__
from cpgmodule.utils import printlog, read_grp_file1


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
            "Data matrix containing beta values. The first row contains "
            "unique sample IDs and the first column contains unique CpG/probe "
            "IDs. Non-numeric values are treated as missing. Plain text, "
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
            "Two groups use Mann-Whitney U; three or more groups use "
            "Kruskal-Wallis."
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


def mann_whitney_test(a, b):
    """
    Perform a two-sided Mann-Whitney U test.

    Returns
    -------
    tuple[float, float]
        p-value and U statistic. NaN values are returned if the test cannot
        be computed.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]

    if a.size == 0 or b.size == 0:
        return math.nan, math.nan

    try:
        result = stats.mannwhitneyu(
            a,
            b,
            alternative="two-sided",
        )
    except ValueError:
        return math.nan, math.nan

    return float(result.pvalue), float(result.statistic)


def kruskal_wallis_test(groups):
    """
    Perform the Kruskal-Wallis H-test.

    Parameters
    ----------
    groups : iterable of array-like
        Group-wise observations.

    Returns
    -------
    tuple[float, float]
        p-value and H statistic. NaN values are returned if the test cannot
        be computed.
    """
    cleaned = []

    for values in groups:
        array = np.asarray(values, dtype=float)
        array = array[np.isfinite(array)]
        if array.size > 0:
            cleaned.append(array)

    if len(cleaned) < 2:
        return math.nan, math.nan

    try:
        result = stats.kruskal(*cleaned)
    except ValueError:
        return math.nan, math.nan

    return float(result.pvalue), float(result.statistic)


# Backward-compatible aliases for code that may import the historical names.
mwu_test = mann_whitney_test


def kruskal_test(*args):
    """Backward-compatible wrapper around :func:`kruskal_wallis_test`."""
    return kruskal_wallis_test(args)


def read_groups(group_file):
    """
    Read and validate the sample-to-group mapping.

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

    if len(group_ids) < 2:
        raise ValueError("At least two groups are required")

    if len(group_ids) == 2:
        printlog("Perform Mann-Whitney U test for two groups ...")
    else:
        printlog("Perform Kruskal-Wallis H-test ...")

    return sample_to_group, group_ids, group_to_samples


def parse_beta(value):
    """Parse one beta value, returning NaN for missing/non-numeric input."""
    try:
        beta = float(value)
    except ValueError:
        return math.nan

    if not math.isfinite(beta):
        return math.nan

    return beta


def analyze_probe(
    beta_values,
    sample_ids,
    sample_to_group,
    group_ids,
):
    """
    Run the appropriate nonparametric test for one CpG/probe.

    Returns
    -------
    tuple[float, float]
        p-value and test statistic.
    """
    group_values = {
        group_id: []
        for group_id in group_ids
    }

    for sample_id, value in zip(sample_ids, beta_values):
        group_id = sample_to_group.get(sample_id)

        if group_id is None:
            continue

        beta = parse_beta(value)
        if math.isfinite(beta):
            group_values[group_id].append(beta)

    if len(group_ids) == 2:
        return mann_whitney_test(
            group_values[group_ids[0]],
            group_values[group_ids[1]],
        )

    return kruskal_wallis_test(
        group_values[group_id]
        for group_id in group_ids
    )


def analyze_file(input_file, sample_to_group, group_ids):
    """
    Analyze all CpG/probe rows in the input matrix.

    Returns
    -------
    tuple
        Header line, original records, probe IDs, p-values, and test statistics.
    """
    header = None
    sample_ids = None

    records = []
    probe_ids = []
    p_values = []
    statistics = []

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
        beta_values = fields[1:]

        if len(beta_values) != len(sample_ids):
            print(
                f"Warning: line {line_number} has {len(beta_values)} data "
                f"fields but header has {len(sample_ids)} samples. Missing "
                "values will be padded with NaN; extra values will be ignored.",
                file=sys.stderr,
            )

        if len(beta_values) < len(sample_ids):
            beta_values = beta_values + ["NaN"] * (
                len(sample_ids) - len(beta_values)
            )
        elif len(beta_values) > len(sample_ids):
            beta_values = beta_values[: len(sample_ids)]

        p_value, statistic = analyze_probe(
            beta_values=beta_values,
            sample_ids=sample_ids,
            sample_to_group=sample_to_group,
            group_ids=group_ids,
        )

        records.append(line)
        probe_ids.append(probe_id)
        p_values.append(p_value)
        statistics.append(statistic)

    if header is None:
        raise ValueError("Input data file is empty")

    if not records:
        raise ValueError("Input data file contains no CpG/probe rows")

    return (
        header,
        records,
        probe_ids,
        p_values,
        statistics,
    )


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
        if math.isfinite(p_value) and 0.0 <= p_value <= 1.0:
            valid_indices.append(index)
            valid_pvalues.append(p_value)

    adjusted = [math.nan] * len(p_values)

    if valid_pvalues:
        corrected = padjust.multiple_testing_correction(valid_pvalues)

        for index, q_value in zip(valid_indices, corrected):
            adjusted[index] = float(q_value)

    return adjusted


def format_stat(value):
    """Format a numeric statistic while preserving missing values."""
    if not math.isfinite(value):
        return "NaN"
    return str(value)


def run_nonparametric(input_file, group_file, output_prefix):
    """
    Run nonparametric differential methylation analysis.

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
        p_values,
        _statistics,
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
            header + "\tpval\tadj.pval",
            file=fout,
        )

        for (
            line,
            probe_id,
            p_value,
            adjusted_pvalue,
        ) in zip(
            records,
            probe_ids,
            p_values,
            adjusted_pvalues,
        ):
            print(
                line
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
        run_nonparametric(
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

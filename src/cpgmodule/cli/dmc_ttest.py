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
Perform differential CpG analysis on beta values using parametric tests.

For exactly two groups:
* Student's t-test is used by default.
* Welch's t-test is used with ``--welch``.
* A paired t-test is used with ``--paired``.

For three or more groups, one-way ANOVA is used.

The first row of the input matrix must contain unique sample IDs, and the
first column must contain CpG/probe IDs. Non-numeric values are treated as
missing.

Output
------
Three columns are appended to the original input table:

* delta_beta
* pval
* adj.pval

For two-group comparisons, ``delta_beta`` is defined as:

    mean(group_1) - mean(group_2)

For multi-group ANOVA, ``delta_beta`` is reported as ``n/a``.

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
            "Two groups use a t-test; three or more groups use ANOVA."
        ),
    )

    parser.add_argument(
        "-p",
        "--paired",
        action="store_true",
        default=False,
        help=(
            "Use a paired t-test for two groups. Samples are paired by their "
            "order within each group."
        ),
    )

    parser.add_argument(
        "-w",
        "--welch",
        action="store_true",
        default=False,
        dest="welch_ttest",
        help=(
            "Use Welch's two-sample t-test, which does not assume equal "
            "variances. Ignored when --paired is used."
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


def standard_ttest(a, b, equal_var=True, nan_policy="omit"):
    """
    Perform an independent two-sample t-test.

    Returns
    -------
    tuple[float, float]
        p-value and t statistic. NaNs are returned if the test cannot be
        computed.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if nan_policy == "omit":
        a = a[np.isfinite(a)]
        b = b[np.isfinite(b)]

    if a.size < 2 or b.size < 2:
        return math.nan, math.nan

    try:
        result = stats.ttest_ind(
            a,
            b,
            equal_var=equal_var,
            nan_policy=nan_policy,
        )
    except ValueError:
        return math.nan, math.nan

    return float(result.pvalue), float(result.statistic)


def paired_ttest(a, b, nan_policy="omit"):
    """
    Perform a paired t-test.

    Missing values are removed pairwise so pairing is preserved.

    Returns
    -------
    tuple[float, float]
        p-value and t statistic. NaNs are returned if the test cannot be
        computed.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    if a.size != b.size:
        return math.nan, math.nan

    if nan_policy == "omit":
        keep = np.isfinite(a) & np.isfinite(b)
        a = a[keep]
        b = b[keep]

    if a.size < 2:
        return math.nan, math.nan

    try:
        result = stats.ttest_rel(
            a,
            b,
            nan_policy=nan_policy,
        )
    except ValueError:
        return math.nan, math.nan

    return float(result.pvalue), float(result.statistic)


def anova(*groups):
    """
    Perform one-way ANOVA across independent groups.

    Returns
    -------
    tuple[float, float]
        p-value and F statistic. NaNs are returned if the test cannot be
        computed.
    """
    cleaned = []

    for values in groups:
        array = np.asarray(values, dtype=float)
        array = array[np.isfinite(array)]

        if array.size == 0:
            return math.nan, math.nan

        cleaned.append(array)

    if len(cleaned) < 2:
        return math.nan, math.nan

    try:
        result = stats.f_oneway(*cleaned)
    except ValueError:
        return math.nan, math.nan

    return float(result.pvalue), float(result.statistic)


def read_groups(group_file, paired=False):
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

    if paired:
        if len(group_ids) != 2:
            raise ValueError(
                "--paired can only be used when exactly two groups are defined"
            )

        if (
            len(group_to_samples[group_ids[0]])
            != len(group_to_samples[group_ids[1]])
        ):
            raise ValueError(
                "Paired t-test requires equal sample counts in the two groups"
            )

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


def build_group_values(
    sample_ids,
    beta_values,
    sample_to_group,
    group_ids,
):
    """Build group-wise beta-value arrays in input-column order."""
    group_values = {
        group_id: []
        for group_id in group_ids
    }

    for sample_id, raw_value in zip(sample_ids, beta_values):
        group_id = sample_to_group.get(sample_id)

        if group_id is None:
            continue

        group_values[group_id].append(parse_beta(raw_value))

    return group_values


def analyze_probe(
    sample_ids,
    beta_values,
    sample_to_group,
    group_ids,
    paired=False,
    welch=False,
):
    """
    Run the appropriate statistical test for one CpG/probe.

    Returns
    -------
    tuple[float, float, float]
        delta_beta, p-value, test statistic.
    """
    group_values = build_group_values(
        sample_ids=sample_ids,
        beta_values=beta_values,
        sample_to_group=sample_to_group,
        group_ids=group_ids,
    )

    if len(group_ids) == 2:
        group1 = np.asarray(
            group_values[group_ids[0]],
            dtype=float,
        )
        group2 = np.asarray(
            group_values[group_ids[1]],
            dtype=float,
        )

        finite1 = group1[np.isfinite(group1)]
        finite2 = group2[np.isfinite(group2)]

        if finite1.size > 0 and finite2.size > 0:
            delta_beta = float(
                np.mean(finite1) - np.mean(finite2)
            )
        else:
            delta_beta = math.nan

        if paired:
            p_value, statistic = paired_ttest(
                group1,
                group2,
                nan_policy="omit",
            )
        else:
            # Welch's t-test means equal_var=False.
            p_value, statistic = standard_ttest(
                group1,
                group2,
                equal_var=not welch,
                nan_policy="omit",
            )

        return delta_beta, p_value, statistic

    groups = [
        np.asarray(group_values[group_id], dtype=float)
        for group_id in group_ids
    ]

    p_value, statistic = anova(*groups)

    return math.nan, p_value, statistic


def analyze_file(
    input_file,
    sample_to_group,
    group_ids,
    paired=False,
    welch=False,
):
    """
    Analyze all CpG/probe rows in the input matrix.

    Returns
    -------
    tuple
        Header, original data rows, probe IDs, delta-beta values, p-values,
        and test statistics.
    """
    header = None
    sample_ids = None

    records = []
    probe_ids = []
    delta_betas = []
    p_values = []
    statistics = []

    warned_missing_samples = False

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
                printlog(
                    "Sample(s) from group file not found in data matrix and "
                    "will be excluded: "
                    + ", ".join(missing)
                )
                warned_missing_samples = True

            if paired:
                present_counts = {
                    group_id: sum(
                        1
                        for sample in sample_ids
                        if sample_to_group.get(sample) == group_id
                    )
                    for group_id in group_ids
                }

                if present_counts[group_ids[0]] != present_counts[group_ids[1]]:
                    raise ValueError(
                        "After matching the data matrix, paired t-test groups "
                        "do not contain the same number of samples"
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

        delta_beta, p_value, statistic = analyze_probe(
            sample_ids=sample_ids,
            beta_values=beta_values,
            sample_to_group=sample_to_group,
            group_ids=group_ids,
            paired=paired,
            welch=welch,
        )

        records.append(line)
        probe_ids.append(probe_id)
        delta_betas.append(delta_beta)
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
        delta_betas,
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
        corrected = padjust.multiple_testing_correction(
            valid_pvalues
        )

        for index, q_value in zip(valid_indices, corrected):
            adjusted[index] = float(q_value)

    return adjusted


def format_delta_beta(value):
    """Format delta-beta output."""
    if not math.isfinite(value):
        return "n/a"
    return str(value)


def format_pvalue(value):
    """Format p-value output."""
    if not math.isfinite(value):
        return "n/a"
    return str(value)


def run_ttest(
    input_file,
    group_file,
    output_prefix,
    paired=False,
    welch=False,
):
    """
    Run t-test/ANOVA differential methylation analysis.

    Returns
    -------
    str
        Path to the generated output file.
    """
    sample_to_group, group_ids, group_to_samples = read_groups(
        group_file,
        paired=paired,
    )

    if len(group_ids) == 2 and paired:
        printlog("Perform paired t-test of two related samples ...")
    elif len(group_ids) == 2 and welch:
        printlog(
            "Perform Welch's t-test of two independent samples ..."
        )
    elif len(group_ids) == 2:
        printlog(
            "Perform Student's t-test of two independent samples ..."
        )
    else:
        printlog("Perform one-way ANOVA ...")

    (
        header,
        records,
        probe_ids,
        delta_betas,
        p_values,
        _statistics,
    ) = analyze_file(
        input_file=input_file,
        sample_to_group=sample_to_group,
        group_ids=group_ids,
        paired=paired,
        welch=welch,
    )

    printlog("Perform Benjamini-Hochberg (FDR) correction ...")
    adjusted_pvalues = adjust_pvalues(p_values)

    output_path = Path(f"{output_prefix}.pval.txt")

    if output_path.parent != Path("."):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    printlog(f'Writing to "{output_path}"')

    with open(output_path, "w") as fout:
        print(
            header + "\tdelta_beta\tpval\tadj.pval",
            file=fout,
        )

        for (
            line,
            probe_id,
            delta_beta,
            p_value,
            adjusted_pvalue,
        ) in zip(
            records,
            probe_ids,
            delta_betas,
            p_values,
            adjusted_pvalues,
        ):
            print(
                line
                + "\t"
                + format_delta_beta(delta_beta)
                + "\t"
                + format_pvalue(p_value)
                + "\t"
                + format_pvalue(adjusted_pvalue),
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

    if args.paired and args.welch_ttest:
        print(
            "Warning: --welch is ignored when --paired is used.",
            file=sys.stderr,
        )

    try:
        run_ttest(
            input_file=args.input_file,
            group_file=args.group_file,
            output_prefix=args.out_file,
            paired=args.paired,
            welch=args.welch_ttest and not args.paired,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

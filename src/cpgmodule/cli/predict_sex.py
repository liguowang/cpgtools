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
Predict biological sex from X-chromosome methylation using the
semi-methylation (SM) ratio.

The method uses the observation that X-chromosome inactivation produces a
higher proportion of semi-methylated CpGs in samples with two X chromosomes.

For each sample, X-linked CpGs are divided into three beta-value ranges:

* low:  0.0 <= beta <= 0.2
* mid:  0.3 <= beta <= 0.7
* high: 0.8 <= beta <= 1.0

The score is:

    log2_SM_ratio = log2(mid / (low + high))

A score greater than the selected cutoff is predicted as Female; a score
lower than the cutoff is predicted as Male. A score exactly equal to the
cutoff, or an undefined score, is reported as Unknown.

Example input
-------------
    CpG_ID    Sample_01    Sample_02    Sample_03
    cg_001    0.831035     0.878022     0.794427
    cg_002    0.249544     0.209949     0.234294
    cg_003    0.845065     0.843957     0.840184

Example output
--------------
    Sample_ID    log2_SM_ratio    Predicted_sex
    Sample_01    -2.249628        Male
    Sample_02    -2.267173        Male
    Sample_03     1.453058        Female
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog


LOW_RANGE = (0.0, 0.2)
MID_RANGE = (0.3, 0.7)
HIGH_RANGE = (0.8, 1.0)


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
            "Tab-separated beta-value matrix. The first row contains sample "
            "IDs and the first column contains CpG IDs."
        ),
    )

    parser.add_argument(
        "-x",
        "--xprobe",
        dest="xprobe_file",
        required=True,
        help=(
            "File containing CpG IDs mapped to chromosome X, one ID per line. "
            "Blank lines and lines beginning with '#' are ignored."
        ),
    )

    parser.add_argument(
        "-c",
        "--cut",
        dest="cutoff",
        type=float,
        default=0.0,
        help=(
            "Cutoff for log2(SM ratio). Values greater than the cutoff are "
            "predicted Female; values below it are predicted Male "
            "[default: %(default)s]."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help=(
            "Output prefix. Predictions are written to "
            "<prefix>.predicted_sex.tsv."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def read_x_probes(xprobe_file):
    """
    Read X-chromosome CpG IDs.

    Returns
    -------
    set[str]
        Unique X-linked CpG IDs.
    """
    printlog(f'Reading X probes from: "{xprobe_file}"')

    probes = set()

    for raw_line in ireader.reader(xprobe_file):
        line = raw_line.strip()

        if not line or line.startswith("#"):
            continue

        probes.add(line)

    if not probes:
        raise ValueError(
            f'No X-chromosome CpG IDs were found in "{xprobe_file}"'
        )

    printlog(f"Total {len(probes)} X probes loaded.")

    return probes


def read_beta_matrix(input_file):
    """
    Read a beta-value matrix.

    Non-numeric values are converted to NaN rather than causing entire CpG
    rows to be discarded.

    Returns
    -------
    pandas.DataFrame
        Numeric beta-value matrix indexed by CpG ID.
    """
    printlog(f'Reading input file: "{input_file}"')

    try:
        frame = pd.read_csv(
            input_file,
            index_col=0,
            sep="\t",
        )
    except Exception as exc:
        raise ValueError(
            f'Cannot read beta-value matrix "{input_file}": {exc}'
        ) from exc

    if frame.empty:
        raise ValueError("Input beta-value matrix is empty")

    if frame.columns.empty:
        raise ValueError("Input beta-value matrix contains no samples")

    if not frame.index.is_unique:
        raise ValueError("CpG IDs in the input matrix are not unique")

    if not frame.columns.is_unique:
        raise ValueError("Sample IDs in the input matrix are not unique")

    numeric = frame.apply(
        pd.to_numeric,
        errors="coerce",
    )

    invalid_count = int(numeric.isna().sum().sum())
    if invalid_count:
        printlog(
            f"{invalid_count} missing/non-numeric beta value(s) will be "
            "ignored on a per-sample basis."
        )

    return numeric


def select_x_probes(beta_matrix, x_probes, input_file=None):
    """
    Select X-linked CpGs found in the beta matrix.

    Returns
    -------
    pandas.DataFrame
        X-linked subset of the beta matrix.
    """
    found = beta_matrix.index.intersection(sorted(x_probes))

    source = f' from file "{input_file}"' if input_file else ""
    printlog(
        f"Found {len(found)} X-linked CpGs{source}."
    )

    if len(found) == 0:
        raise ValueError(
            "None of the supplied X-chromosome CpG IDs were found in the "
            "input beta-value matrix"
        )

    return beta_matrix.loc[found]


def _count_range(values, lower, upper):
    """Count finite beta values within an inclusive interval."""
    array = np.asarray(values, dtype=float)
    valid = np.isfinite(array)

    return int(
        np.sum(
            valid
            & (array >= lower)
            & (array <= upper)
        )
    )


def calculate_sm_ratio(
    values,
    low_range=LOW_RANGE,
    mid_range=MID_RANGE,
    high_range=HIGH_RANGE,
):
    """
    Calculate log2 semi-methylation ratio for one sample.

    Returns
    -------
    tuple
        log2 ratio, low count, mid count, high count, valid X-CpG count.
    """
    array = np.asarray(values, dtype=float)
    valid_count = int(np.sum(np.isfinite(array)))

    low_count = _count_range(
        array,
        low_range[0],
        low_range[1],
    )
    mid_count = _count_range(
        array,
        mid_range[0],
        mid_range[1],
    )
    high_count = _count_range(
        array,
        high_range[0],
        high_range[1],
    )

    denominator = low_count + high_count

    # The historical method is undefined when either the numerator or
    # denominator is zero. Report NaN instead of relying on divide-by-zero
    # warnings or infinities.
    if mid_count == 0 or denominator == 0:
        ratio = math.nan
    else:
        ratio = math.log2(mid_count / denominator)

    return (
        ratio,
        low_count,
        mid_count,
        high_count,
        valid_count,
    )


def predict_from_ratio(ratio, cutoff=0.0):
    """Convert a log2 SM ratio to a sex prediction."""
    if not math.isfinite(ratio):
        return "Unknown"

    if ratio > cutoff:
        return "Female"

    if ratio < cutoff:
        return "Male"

    return "Unknown"


def predict_samples(beta_matrix, cutoff=0.0):
    """
    Predict sex for each sample in an X-linked beta matrix.

    Returns
    -------
    pandas.DataFrame
        Per-sample SM-ratio calculations and predictions.
    """
    rows = []

    for sample_id in beta_matrix.columns:
        (
            ratio,
            low_count,
            mid_count,
            high_count,
            valid_count,
        ) = calculate_sm_ratio(beta_matrix[sample_id].to_numpy())

        rows.append(
            {
                "Sample_ID": sample_id,
                "log2_SM_ratio": ratio,
                "Predicted_sex": predict_from_ratio(
                    ratio,
                    cutoff=cutoff,
                ),
                "X_CpGs_used": valid_count,
                "Low_beta_CpGs": low_count,
                "Mid_beta_CpGs": mid_count,
                "High_beta_CpGs": high_count,
            }
        )

    return pd.DataFrame(rows).set_index("Sample_ID")


def run_predict_sex(
    input_file,
    xprobe_file,
    output_prefix,
    cutoff=0.0,
):
    """
    Run X-chromosome semi-methylation sex prediction.

    Returns
    -------
    str
        Path to the generated prediction table.
    """
    x_probes = read_x_probes(xprobe_file)

    beta_matrix = read_beta_matrix(input_file)

    x_beta_matrix = select_x_probes(
        beta_matrix,
        x_probes,
        input_file=input_file,
    )

    predictions = predict_samples(
        x_beta_matrix,
        cutoff=cutoff,
    )

    output_path = Path(
        f"{output_prefix}.predicted_sex.tsv"
    )

    if output_path.parent != Path("."):
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    printlog(f'Writing to file: "{output_path}"')

    predictions.to_csv(
        output_path,
        sep="\t",
        index_label="Sample_ID",
        na_rep="NaN",
    )

    return str(output_path)


def main(argv=None):
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str] or None
        Optional command-line argument list. When None, argparse reads
        ``sys.argv``. Passing a list allows this function to be called directly
        from the CpGtools dispatcher or from tests.

    Returns
    -------
    int
        Process-style return code.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        run_predict_sex(
            input_file=args.input_file,
            xprobe_file=args.xprobe_file,
            output_prefix=args.out_file,
            cutoff=args.cutoff,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

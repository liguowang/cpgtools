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
Perform differential CpG analysis using logistic regression on methylated and
total read counts.

Input methylation values are represented as:

    methylated_count,total_count

where both values are non-negative integers and
``methylated_count <= total_count``.

The first row of the data matrix contains sample IDs and the first column
contains CpG/probe IDs. The group/covariate file is parsed with
``cpgmodule.utils.read_grp_file2``. The first covariate is treated as the
primary grouping variable, while additional variables are included as
covariates.

The statistical model is fitted in R using:

    glm(cbind(m, t - m) ~ covariates, family = ...)

Supported families
------------------
1. quasibinomial (default)
2. binomial

Output files
------------
<prefix>.r
    Generated R script.

<prefix>.results.txt
    Full coefficient and p-value table written by R.

<prefix>.warnings.txt
    R stderr output.

<prefix>.pval.txt
    Original input table with primary-variable p-value and adjusted p-value.
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

from cpgmodule import ireader
from cpgmodule import padjust
from cpgmodule._version import __version__
from cpgmodule.utils import printlog, read_grp_file2


COUNT_RE = re.compile(r"^(\d+)\s*,\s*(\d+)$")
FAMILIES = {
    1: "quasibinomial",
    2: "binomial",
}


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
            "Data matrix containing methylated,total counts. The first row "
            "contains unique sample IDs and the first column contains unique "
            "CpG/probe IDs. Plain text, .gz, .bz2, and other formats "
            "supported by cpgmodule.ireader may be used."
        ),
    )

    parser.add_argument(
        "-g",
        "--group",
        dest="group_file",
        required=True,
        help=(
            "Group/covariate file. The first variable is the primary grouping "
            "variable and must be categorical. Additional variables may be "
            "categorical or continuous."
        ),
    )

    parser.add_argument(
        "-f",
        "--family",
        dest="family_func",
        type=int,
        choices=sorted(FAMILIES),
        default=1,
        help=(
            "GLM family: 1=quasibinomial, 2=binomial "
            "[default: %(default)s]."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help=(
            "Output prefix. Files include <prefix>.r, "
            "<prefix>.results.txt, <prefix>.warnings.txt, and "
            "<prefix>.pval.txt."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def r_quote(value):
    """Return a safely quoted R string literal."""
    value = str(value)
    value = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{value}"'


def parse_count(value):
    """
    Parse one methylated_count,total_count field.

    Returns
    -------
    tuple[int, int] or None
        Parsed counts, or None for missing/non-count text.

    Raises
    ------
    ValueError
        If counts are syntactically valid but biologically impossible.
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


def validate_group_data(samples, covariate_names, covariates, variable_types):
    """Validate metadata returned by read_grp_file2()."""
    if not samples:
        raise ValueError("Group file contains no samples")

    if not covariate_names:
        raise ValueError(
            "Group file must contain at least one grouping variable"
        )

    primary = covariate_names[0]
    if variable_types.get(primary) != "categorical":
        raise ValueError(
            f'Primary grouping variable "{primary}" must be categorical'
        )

    for name in covariate_names:
        if name not in variable_types:
            raise ValueError(
                f'Cannot determine variable type for "{name}"'
            )

        missing = [
            sample
            for sample in samples
            if sample not in covariates[name]
        ]
        if missing:
            raise ValueError(
                f'Covariate "{name}" is missing values for: '
                + ", ".join(missing)
            )


def write_r_functions(
    rout,
    results_file,
    covariate_names,
    family_name,
):
    """Write reusable logistic-regression helper code to the R script."""
    cov_args = ",".join(covariate_names)
    formula = "+".join(covariate_names)
    results_literal = r_quote(results_file)

    print(
        f"logit_fit <- function(cgid, m, t, {cov_args}, append=TRUE) {{",
        file=rout,
    )
    print(
        "\tdat <- data.frame("
        "m=m, t=t, "
        + ", ".join(f"{name}={name}" for name in covariate_names)
        + ")",
        file=rout,
    )
    print(
        f"\tfit <- try(glm(cbind(m, t - m) ~ {formula}, "
        f"family={family_name}, data=na.omit(dat)), silent=TRUE)",
        file=rout,
    )
    print("\tif (inherits(fit, \"try-error\")) {", file=rout)
    print("\t\treturn(invisible(NULL))", file=rout)
    print("\t}", file=rout)

    print("\ttab <- coef(summary(fit))", file=rout)
    print("\tcoefs <- tab[,1]", file=rout)
    print("\tpvals <- tab[,4]", file=rout)
    print(
        "\tif (any(pvals > 1, na.rm=TRUE)) { pvals <- pvals + NA }",
        file=rout,
    )
    print(
        "\tif (sum(m, na.rm=TRUE) == 0) { pvals <- pvals + NA }",
        file=rout,
    )

    print("\tnms <- row.names(tab)", file=rout)
    print('\tnms <- gsub("2", "", nms)', file=rout)

    print(
        "\trow <- matrix(c(cgid, as.vector(coefs), "
        "as.vector(pvals)), nrow=1)",
        file=rout,
    )

    print("\tif (!append) {", file=rout)
    print(
        f"\t\twrite.table(file={results_literal}, x=row, quote=FALSE, "
        "row.names=FALSE, sep=\"\\t\", "
        'col.names=c("ID", paste(nms, "coef", sep="."), '
        'paste(nms, "pval", sep=".")))',
        file=rout,
    )
    print("\t} else {", file=rout)
    print(
        f"\t\twrite.table(file={results_literal}, x=row, quote=FALSE, "
        "row.names=FALSE, sep=\"\\t\", col.names=FALSE, append=TRUE)",
        file=rout,
    )
    print("\t}", file=rout)
    print("}", file=rout)
    print("", file=rout)


def write_covariates(
    rout,
    sample_ids,
    covariate_names,
    covariates,
    variable_types,
):
    """Write R covariate vectors in the data-matrix sample order."""
    for name in covariate_names:
        values = [
            covariates[name].get(sample)
            for sample in sample_ids
        ]

        if variable_types[name] == "continuous":
            encoded = ",".join(
                "NA" if value is None else str(value)
                for value in values
            )
            print(f"{name} <- c({encoded})", file=rout)

        elif variable_types[name] == "categorical":
            encoded = ",".join(
                "NA" if value is None else r_quote(value)
                for value in values
            )
            print(
                f"{name} <- as.factor(c({encoded}))",
                file=rout,
            )

        else:
            raise ValueError(
                f'Unknown variable type for "{name}"'
            )

    print("", file=rout)


def write_analysis_calls(
    rout,
    input_file,
    samples,
    covariate_names,
    covariates,
    variable_types,
):
    """
    Parse the count matrix and write one R GLM call per CpG/probe.

    Returns
    -------
    tuple[int, str]
        Number of probes written and original header line.
    """
    sample_ids = None
    header_line = None
    probe_count = 0

    sample_set = set(samples)

    for line_number, raw_line in enumerate(
        ireader.reader(input_file),
        start=1,
    ):
        line = raw_line.rstrip("\r\n")

        if not line:
            continue

        fields = line.split()

        if sample_ids is None:
            if len(fields) < 2:
                raise ValueError(
                    "Input matrix header must contain at least one sample ID"
                )

            header_line = line
            sample_ids = fields[1:]

            if len(sample_ids) != len(set(sample_ids)):
                raise ValueError(
                    "Sample IDs in input matrix are not unique"
                )

            missing = [
                sample
                for sample in samples
                if sample not in sample_ids
            ]
            if missing:
                raise ValueError(
                    "Cannot find sample ID(s) in data file: "
                    + ", ".join(missing)
                )

            unknown = [
                sample
                for sample in sample_ids
                if sample not in sample_set
            ]
            if unknown:
                printlog(
                    "Ignoring sample(s) not present in group file: "
                    + ", ".join(unknown)
                )

            write_covariates(
                rout,
                sample_ids,
                covariate_names,
                covariates,
                variable_types,
            )
            continue

        probe_id = fields[0]
        count_fields = fields[1:]

        if len(count_fields) != len(sample_ids):
            print(
                f"Warning: line {line_number} has {len(count_fields)} data "
                f"fields but header has {len(sample_ids)} samples. Missing "
                "values will be padded with NA; extra values will be ignored.",
                file=sys.stderr,
            )

        methylated_reads = []
        total_reads = []

        for sample_id, value in zip(
            sample_ids,
            count_fields[: len(sample_ids)],
        ):
            if sample_id not in sample_set:
                methylated_reads.append("NA")
                total_reads.append("NA")
                continue

            try:
                parsed = parse_count(value)
            except ValueError as exc:
                raise ValueError(
                    f'Invalid count data on line {line_number}, '
                    f'probe "{probe_id}", sample "{sample_id}": {exc}'
                ) from exc

            if parsed is None:
                methylated_reads.append("NA")
                total_reads.append("NA")
            else:
                methylated, total = parsed
                methylated_reads.append(str(methylated))
                total_reads.append(str(total))

        while len(methylated_reads) < len(sample_ids):
            methylated_reads.append("NA")
            total_reads.append("NA")

        probe_count += 1
        append_flag = "FALSE" if probe_count == 1 else "TRUE"

        print(
            "logit_fit("
            f"{r_quote(probe_id)}, "
            f"c({','.join(methylated_reads)}), "
            f"c({','.join(total_reads)}), "
            f"{','.join(covariate_names)}, "
            f"append={append_flag}"
            ")",
            file=rout,
        )

    if sample_ids is None:
        raise ValueError("Input data file is empty")

    if probe_count == 0:
        raise ValueError("Input data file contains no CpG/probe rows")

    return probe_count, header_line


def generate_r_script(
    input_file,
    group_file,
    output_prefix,
    family_func=1,
):
    """
    Generate the R script for CpG-wise logistic regression.

    Returns
    -------
    tuple
        R script path, results path, warnings path, primary variable,
        probe count, and original header.
    """
    family_name = FAMILIES[family_func]

    printlog(f'Read group file "{group_file}" ...')

    samples, covariate_names, covariates, variable_types = read_grp_file2(
        group_file
    )

    validate_group_data(
        samples,
        covariate_names,
        covariates,
        variable_types,
    )

    for name in covariate_names:
        printlog(f"{name}: {variable_types[name]}")

    primary_variable = covariate_names[0]

    output_prefix = Path(output_prefix)
    if output_prefix.parent != Path("."):
        output_prefix.parent.mkdir(parents=True, exist_ok=True)

    r_script = f"{output_prefix}.r"
    results_file = f"{output_prefix}.results.txt"
    warnings_file = f"{output_prefix}.warnings.txt"

    printlog(f'Processing file "{input_file}" ...')

    with open(r_script, "w") as rout:
        write_r_functions(
            rout,
            results_file,
            covariate_names,
            family_name,
        )

        probe_count, header_line = write_analysis_calls(
            rout,
            input_file,
            samples,
            covariate_names,
            covariates,
            variable_types,
        )

    return (
        r_script,
        results_file,
        warnings_file,
        primary_variable,
        probe_count,
        header_line,
    )


def run_rscript(r_script, warnings_file):
    """Execute the generated R script and capture stderr."""
    rscript = shutil.which("Rscript")

    if rscript is None:
        raise RuntimeError(
            "Rscript was not found in PATH. Install R and ensure Rscript "
            "is available before running dmc_logit."
        )

    printlog(f'Running R script "{r_script}" ...')

    with open(warnings_file, "w") as err:
        result = subprocess.run(
            [rscript, r_script],
            stderr=err,
            text=True,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"Rscript failed with exit status {result.returncode}. "
            f'See "{warnings_file}" for details.'
        )


def find_primary_pvalue_column(header, primary_variable):
    """
    Find the p-value column corresponding to the primary grouping variable.
    """
    candidates = [
        index
        for index, name in enumerate(header)
        if name.startswith(primary_variable) and name.endswith(".pval")
    ]

    if not candidates:
        raise ValueError(
            f'Cannot find p-value column for primary variable '
            f'"{primary_variable}" in R results'
        )

    return candidates[0]


def read_primary_pvalues(results_file, primary_variable):
    """
    Read primary-variable p-values from the R results table.

    Returns
    -------
    dict[str, float]
        Probe ID to valid p-value mapping.
    """
    results_path = Path(results_file)

    if not results_path.exists() or results_path.stat().st_size == 0:
        raise RuntimeError(
            f'R did not produce a valid results file: "{results_file}"'
        )

    pvalues = {}

    with open(results_path, "r", newline="") as fin:
        reader = csv.reader(fin, delimiter="\t")

        try:
            header = next(reader)
        except StopIteration as exc:
            raise RuntimeError(
                f'R results file is empty: "{results_file}"'
            ) from exc

        pvalue_index = find_primary_pvalue_column(
            header,
            primary_variable,
        )

        for row in reader:
            if not row:
                continue

            probe_id = row[0]

            if pvalue_index >= len(row):
                continue

            try:
                pvalue = float(row[pvalue_index])
            except ValueError:
                continue

            if 0.0 <= pvalue <= 1.0:
                pvalues[probe_id] = pvalue

    return pvalues


def adjust_primary_pvalues(pvalues):
    """
    Apply CpGtools multiple-testing correction.

    Returns
    -------
    dict[str, tuple[float, float]]
        Probe ID mapped to (p-value, adjusted p-value).
    """
    probe_ids = list(pvalues)
    raw = [pvalues[probe_id] for probe_id in probe_ids]

    if not raw:
        return {}

    adjusted = padjust.multiple_testing_correction(raw)

    return {
        probe_id: (float(pvalue), float(qvalue))
        for probe_id, pvalue, qvalue in zip(
            probe_ids,
            raw,
            adjusted,
        )
    }


def write_pvalue_table(input_file, output_file, adjusted_pvalues):
    """Append primary p-value and adjusted p-value to the input table."""
    printlog(f'Writing to "{output_file}"')

    wrote_header = False

    with open(output_file, "w") as fout:
        for raw_line in ireader.reader(input_file):
            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            if not wrote_header:
                print(
                    line + "\tpval\tadj.pval",
                    file=fout,
                )
                wrote_header = True
                continue

            fields = line.split()
            if not fields:
                continue

            probe_id = fields[0]

            if probe_id in adjusted_pvalues:
                pvalue, qvalue = adjusted_pvalues[probe_id]
                print(
                    f"{line}\t{pvalue}\t{qvalue}",
                    file=fout,
                )
            else:
                print(
                    f"{line}\tNaN\tNaN",
                    file=fout,
                )


def run_logit(
    input_file,
    group_file,
    output_prefix,
    family_func=1,
):
    """
    Run logistic-regression differential methylation analysis.

    Returns
    -------
    str
        Path to the final p-value table.
    """
    (
        r_script,
        results_file,
        warnings_file,
        primary_variable,
        probe_count,
        _,
    ) = generate_r_script(
        input_file=input_file,
        group_file=group_file,
        output_prefix=output_prefix,
        family_func=family_func,
    )

    printlog(f"Generated R commands for {probe_count} CpG/probe(s)")

    run_rscript(
        r_script=r_script,
        warnings_file=warnings_file,
    )

    printlog("Perform Benjamini-Hochberg (FDR) correction ...")

    pvalues = read_primary_pvalues(
        results_file=results_file,
        primary_variable=primary_variable,
    )

    adjusted = adjust_primary_pvalues(pvalues)

    output_file = f"{output_prefix}.pval.txt"

    write_pvalue_table(
        input_file=input_file,
        output_file=output_file,
        adjusted_pvalues=adjusted,
    )

    return output_file


def main(argv=None):
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str] or None
        Optional command-line argument list. When None, argparse reads
        ``sys.argv``. Passing a list allows the function to be called directly
        from the CpGtools dispatcher or from tests.

    Returns
    -------
    int
        Process-style return code.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        run_logit(
            input_file=args.input_file,
            group_file=args.group_file,
            output_prefix=args.out_file,
            family_func=args.family_func,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

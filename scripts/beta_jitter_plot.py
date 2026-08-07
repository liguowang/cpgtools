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

"""Generate jitter and bean plots for DNA methylation beta values.

The input must be a tab-delimited matrix whose first column contains CpG IDs and
whose remaining columns contain sample beta values. A random subset of CpGs may
be used to reduce plot density. Plotting is performed with R and the ``beanplot``
package.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, Sequence

import pandas as pd

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


class JitterPlotError(RuntimeError):
    """Raised when the input data or plotting process is invalid."""


def r_quote(value: str) -> str:
    """Return a double-quoted R string with backslashes and quotes escaped."""
    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def read_beta_matrix(path: Path) -> pd.DataFrame:
    """Read and validate a DNA methylation beta-value matrix."""
    if not path.is_file():
        raise JitterPlotError(f"Input file does not exist: {path}")

    try:
        data = pd.read_csv(path, sep="\t", compression="infer")
    except Exception as exc:  # pragma: no cover - pandas provides details
        raise JitterPlotError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] < 2:
        raise JitterPlotError(
            "The input must contain a CpG-ID column and at least one sample column."
        )
    if data.empty:
        raise JitterPlotError("The input matrix contains no CpG rows.")

    sample_names = data.columns[1:].astype(str)
    if sample_names.duplicated().any():
        duplicates = sorted(set(sample_names[sample_names.duplicated()]))
        raise JitterPlotError(
            "Duplicate sample IDs were found: " + ", ".join(duplicates)
        )

    numeric = data.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    if numeric.notna().sum().sum() == 0:
        raise JitterPlotError("No numeric beta values were found in the sample columns.")

    invalid = ((numeric < 0.0) | (numeric > 1.0)).sum().sum()
    if invalid:
        printlog(
            f"Warning: {invalid} beta values fall outside the expected [0, 1] range."
        )

    result = pd.concat([data.iloc[:, [0]], numeric], axis=1)
    result.columns = [str(column) for column in data.columns]
    return result


def select_plot_data(data: pd.DataFrame, fraction: float, seed: int) -> pd.DataFrame:
    """Select the CpG rows used for plotting."""
    if fraction >= 1.0:
        return data

    number_rows = max(1, round(len(data) * fraction))
    return data.sample(n=number_rows, replace=False, random_state=seed)


def write_r_script(
    script_path: Path,
    data_path: Path,
    image_path: Path,
    sample_names: Sequence[str],
    width: int,
    height: int,
    point_size: float,
    jitter: float,
) -> None:
    """Write the R script used to generate the jitter and bean plot."""
    sample_vector = ", ".join(r_quote(name) for name in sample_names)

    script = f"""suppressPackageStartupMessages(library(beanplot))

d <- read.delim(
    file = {r_quote(str(data_path.resolve()))},
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sample_names <- c({sample_vector})
missing_samples <- setdiff(sample_names, colnames(d))
if (length(missing_samples) > 0) {{
    stop(paste("Samples missing from plot-data file:", paste(missing_samples, collapse=", ")))
}}

values <- lapply(sample_names, function(name) {{
    x <- suppressWarnings(as.numeric(d[[name]]))
    x[is.finite(x)]
}})
names(values) <- sample_names

if (any(lengths(values) == 0)) {{
    empty_samples <- names(values)[lengths(values) == 0]
    stop(paste("No finite beta values for sample(s):", paste(empty_samples, collapse=", ")))
}}

png(
    filename = {r_quote(str(image_path.resolve()))},
    width = {width},
    height = {height},
    units = "px"
)

par(mar = c(8, 4, 3, 1) + 0.1)
stripchart(
    values,
    cex = {point_size},
    col = "#abd9e9",
    vertical = TRUE,
    method = "jitter",
    ylab = "Beta value",
    las = 2,
    jitter = {jitter},
    cex.names = 0.8,
    ylim = c(0, 1)
)
beanplot(
    values,
    cutmin = 0,
    cutmax = 1,
    border = "#d01c8b",
    what = c(1, 1, 1, 0),
    col = NA,
    las = 2,
    add = TRUE
)

dev.off()
"""
    script_path.write_text(script, encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    """Construct the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate a jitter plot with an overlaid bean plot for each sample "
            "in a DNA methylation beta-value matrix."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Tab-delimited beta-value matrix; the first column contains CpG IDs "
            "and later columns contain samples."
        ),
    )
    parser.add_argument(
        "-f",
        "--fraction",
        type=float,
        default=0.5,
        help="Fraction of CpG rows to plot; must be in the interval (0, 1].",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=999,
        help="Random seed used when sampling CpG rows.",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=800,
        help="Output PNG width in pixels.",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=480,
        help="Output PNG height in pixels.",
    )
    parser.add_argument(
        "--point_size",
        type=float,
        default=0.1,
        help="Point-size multiplier for the jitter plot.",
    )
    parser.add_argument(
        "--jitter",
        type=float,
        default=0.3,
        help="Horizontal jitter width used by R's stripchart function.",
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
        "--keep_plot_data",
        action="store_true",
        help="Keep the tab-delimited CpG subset used to generate the plot.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the command-line program."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not 0.0 < args.fraction <= 1.0:
        parser.error("--fraction must be greater than 0 and no greater than 1")
    if args.width < 1 or args.height < 1:
        parser.error("--width and --height must be positive integers")
    if args.point_size <= 0:
        parser.error("--point_size must be greater than 0")
    if args.jitter < 0:
        parser.error("--jitter must be non-negative")

    rscript = shutil.which("Rscript")
    if rscript is None:
        parser.error("Rscript was not found in PATH")

    out_prefix: Path = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    script_path = Path(f"{out_prefix}.r")
    data_path = Path(f"{out_prefix}.plot_data.tsv")
    image_path = Path(f"{out_prefix}.png")

    try:
        printlog(f'Reading beta file: "{args.input_file}"')
        data = read_beta_matrix(args.input_file)
        plot_data = select_plot_data(data, args.fraction, args.seed)

        if len(plot_data) < len(data):
            printlog(
                f"Selected {len(plot_data):,} of {len(data):,} CpGs "
                f"({100.0 * len(plot_data) / len(data):.2f}%) for plotting."
            )
        else:
            printlog(f"Using all {len(plot_data):,} CpGs for plotting.")

        plot_data.to_csv(data_path, sep="\t", index=False, na_rep="NA")
        write_r_script(
            script_path=script_path,
            data_path=data_path,
            image_path=image_path,
            sample_names=list(plot_data.columns[1:]),
            width=args.width,
            height=args.height,
            point_size=args.point_size,
            jitter=args.jitter,
        )

        printlog(f'Running R script: "{script_path}"')
        completed = subprocess.run(
            [rscript, str(script_path)],
            check=False,
            text=True,
            capture_output=True,
        )
        if completed.returncode != 0:
            details = completed.stderr.strip() or completed.stdout.strip()
            raise JitterPlotError(
                f"R plotting failed with exit code {completed.returncode}: {details}"
            )
        if not image_path.is_file():
            raise JitterPlotError(f"R completed but did not create {image_path}")

        printlog(f'Plot saved to: "{image_path}"')
        if not args.keep_plot_data:
            data_path.unlink(missing_ok=True)

    except JitterPlotError as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

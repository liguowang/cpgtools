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

"""Generate a gene-centered CpG density profile.

The program summarizes CpG density across three transcript-oriented regions:

* upstream of the transcription start site (TSS)
* gene body from TSS to transcription end site (TES)
* downstream of the TES

Each region is normalized into bins by ``density_over_range()`` (currently
100 bins per region in CpGtools), producing a 5'-to-3' metagene profile.

The reference gene model is extended using CpGtools' regulatory-domain
functions so that upstream and downstream regions do not extend through
neighboring genes.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd

from cpgmodule import extend_bed
from cpgmodule._version import __version__
from cpgmodule.utils import density_over_range, printlog, read_CpG_bed


class CpGDensityError(RuntimeError):
    """Raised when CpG-density profiling fails."""


REGION_ORDER = ("Upstream", "GeneBody", "Downstream")


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_density_gene_centered",
        description=(
            "Calculate a transcript-oriented CpG density profile across "
            "upstream, gene-body, and downstream regions."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input BED3+ file containing CpG genomic positions. "
            "Compressed input is supported by the CpGtools reader."
        ),
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="gene_file",
        required=True,
        type=Path,
        help="Reference gene model in BED6+ format.",
    )
    parser.add_argument(
        "-u",
        "--upstream",
        dest="upstream_size",
        type=int,
        default=2000,
        help="Maximum upstream extension from the TSS, in base pairs.",
    )
    parser.add_argument(
        "-d",
        "--downstream",
        dest="downstream_size",
        type=int,
        default=2000,
        help="Maximum downstream extension from the TES, in base pairs.",
    )
    parser.add_argument(
        "-c",
        "--min_gene_size",
        "--SizeCut",
        dest="minimum_gene_size",
        type=int,
        default=200,
        help=(
            "Minimum transcript span from TSS to TES. Genes shorter than this "
            "threshold are excluded."
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
        "--format",
        choices=("png", "pdf", "both"),
        default="pdf",
        help="Plot output format.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution of PNG output.",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=10.0,
        help="Plot width in inches.",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=5.0,
        help="Plot height in inches.",
    )
    parser.add_argument(
        "--no_plot",
        action="store_true",
        help="Write the density table without generating a plot.",
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

    if not args.gene_file.is_file():
        parser.error(f"Reference gene file does not exist: {args.gene_file}")

    if args.upstream_size < 0:
        parser.error("--upstream must be non-negative")

    if args.downstream_size < 0:
        parser.error("--downstream must be non-negative")

    if args.minimum_gene_size < 1:
        parser.error("--min_gene_size must be at least 1")

    if args.dpi < 1:
        parser.error("--dpi must be at least 1")

    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def generate_density_profile(
    input_file: Path,
    gene_file: Path,
    upstream_size: int,
    downstream_size: int,
    minimum_gene_size: int,
) -> pd.DataFrame:
    """Calculate CpG density across upstream, gene-body, and downstream bins."""
    printlog(f'Reading CpG file: "{input_file}"')
    cpg_ranges = read_CpG_bed(str(input_file))

    printlog(f'Reading reference gene model: "{gene_file}"')
    basal_ranges = extend_bed.getBasalDomains(str(gene_file))

    extended_ranges = extend_bed.geteExtendedDomains(
        basal_ranges=basal_ranges,
        bedfile=str(gene_file),
        up_ext=upstream_size,
        down_ext=downstream_size,
        min_gene=minimum_gene_size,
        printit=False,
    )

    printlog("Calculating gene-centered CpG density ...")
    up_density, gene_density, down_density = density_over_range(
        extended_ranges,
        cpg_ranges,
    )

    density_sets = {
        "Upstream": up_density,
        "GeneBody": gene_density,
        "Downstream": down_density,
    }

    rows = []
    for region_name in REGION_ORDER:
        values = density_sets[region_name]

        if not values:
            printlog(f'Warning: no density values were produced for "{region_name}".')
            continue

        for position in sorted(values):
            rows.append(
                {
                    "Group": region_name,
                    "Position": int(position),
                    "CpG_count": float(values[position]),
                }
            )

    if not rows:
        raise CpGDensityError(
            "No CpG density values could be generated. Check the CpG and gene files."
        )

    return pd.DataFrame(rows)


def plot_density(
    profile: pd.DataFrame,
    output_path: Path,
    *,
    dpi: int,
    width: float,
    height: float,
) -> None:
    """Plot the concatenated gene-centered CpG density profile."""
    x_values = []
    y_values = []
    tick_positions = []
    boundaries = []
    cursor = 0

    for region_name in REGION_ORDER:
        region = profile.loc[profile["Group"] == region_name].sort_values("Position")
        if region.empty:
            continue

        values = region["CpG_count"].to_numpy(dtype=float)
        start = cursor

        x_values.extend(range(cursor, cursor + len(values)))
        y_values.extend(values.tolist())
        cursor += len(values)

        tick_positions.append((start + cursor - 1) / 2.0)
        boundaries.append(cursor - 0.5)

    if not y_values:
        raise CpGDensityError("No finite CpG density values are available for plotting.")

    figure, axis = plt.subplots(figsize=(width, height))
    axis.plot(x_values, y_values, linewidth=1.5)

    for boundary in boundaries[:-1]:
        axis.axvline(boundary, linestyle="--", linewidth=0.8)

    axis.set_xlim(min(x_values), max(x_values))
    axis.set_ylabel("CpG count")
    axis.set_xlabel("Transcript-oriented genomic regions (5' to 3')")
    axis.set_xticks(tick_positions)
    axis.set_xticklabels(
        [region for region in REGION_ORDER if region in set(profile["Group"])]
    )
    axis.set_title("Gene-centered CpG density profile")

    figure.tight_layout()
    save_kwargs = {"dpi": dpi} if output_path.suffix.lower() == ".png" else {}
    figure.savefig(output_path, **save_kwargs)
    plt.close(figure)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    table_path = Path(f"{args.out_prefix}.CpG_density.tsv")

    try:
        profile = generate_density_profile(
            input_file=args.input_file,
            gene_file=args.gene_file,
            upstream_size=args.upstream_size,
            downstream_size=args.downstream_size,
            minimum_gene_size=args.minimum_gene_size,
        )

        printlog(f'Writing CpG density table: "{table_path}"')
        profile.to_csv(
            table_path,
            sep="\t",
            index=False,
            float_format="%.6f",
        )

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)

            for extension in extensions:
                plot_path = Path(f"{args.out_prefix}.CpG_density.{extension}")
                printlog(f'Writing CpG density plot: "{plot_path}"')
                plot_density(
                    profile=profile,
                    output_path=plot_path,
                    dpi=args.dpi,
                    width=args.width,
                    height=args.height,
                )

    except CpGDensityError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

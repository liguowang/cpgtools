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

"""Generate a gene-centered DNA methylation profile.

The program calculates average beta values across normalized positions in eight
transcript-oriented genomic regions:

* upstream intergenic region
* 5' UTR exons
* coding exons
* first introns
* internal introns
* last introns
* 3' UTR exons
* downstream intergenic region

The input methylation file must be BED6 or BED6+ with beta values in column 5
and strand information in column 6. The reference gene model must be BED12.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Callable, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cpgmodule import BED
from cpgmodule._version import __version__
from cpgmodule.utils import coverage_over_range, printlog, read_CpG_bed


class GeneProfileError(RuntimeError):
    """Raised when input validation or profile generation fails."""


REGION_LABELS = {
    "Upstream_intergenic": "Upstream",
    "Five_prime_UTR": "5' UTR exon",
    "Coding_exon": "Coding exon",
    "First_intron": "First intron",
    "Internal_intron": "Internal intron",
    "Last_intron": "Last intron",
    "Three_prime_UTR": "3' UTR exon",
    "Downstream_intergenic": "Downstream",
}


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_profile_gene_centered",
        description="Generate a transcript-oriented, gene-centered methylation profile.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input BED6+ methylation file. Columns 1-6 are chromosome, start, "
            "end, name, beta value, and strand. Compressed input is supported "
            "when supported by CpGtools' reader."
        ),
    )
    parser.add_argument(
        "-r",
        "--refgene",
        required=True,
        type=Path,
        help="Reference gene model in BED12 format.",
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=2000,
        help="Length of the upstream intergenic region in base pairs.",
    )
    parser.add_argument(
        "-d",
        "--downstream",
        type=int,
        default=2000,
        help="Length of the downstream intergenic region in base pairs.",
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
        help="Write the profile table without generating a plot.",
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
        parser.error(f"Input methylation file does not exist: {args.input_file}")
    if not args.refgene.is_file():
        parser.error(f"Reference gene file does not exist: {args.refgene}")
    if args.upstream < 0:
        parser.error("--upstream must be non-negative")
    if args.downstream < 0:
        parser.error("--downstream must be non-negative")
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def compute_region_profile(
    region_name: str,
    region_ranges,
    cpg_ranges,
) -> pd.DataFrame:
    """Calculate a normalized methylation profile for one genomic region."""
    values = coverage_over_range(region_ranges, cpg_ranges)
    if not values:
        printlog(f'Warning: no methylation values were found for "{region_name}".')
        return pd.DataFrame(
            columns=["Group", "Relative_position(5'->3')", "Average_beta"]
        )

    rows = [
        {
            "Group": region_name,
            "Relative_position(5'->3')": int(position),
            "Average_beta": float(values[position]),
        }
        for position in sorted(values)
    ]
    return pd.DataFrame(rows)


def generate_profiles(
    input_file: Path,
    gene_file: Path,
    upstream_size: int,
    downstream_size: int,
) -> pd.DataFrame:
    """Generate profiles for all supported transcript-oriented regions."""
    printlog(f'Reading CpG file: "{input_file}"')
    cpg_ranges = read_CpG_bed(str(input_file))

    printlog(f'Reading reference gene model: "{gene_file}"')
    ref_gene = BED.ParseBED(str(gene_file))

    region_builders: list[tuple[str, Callable[[], object]]] = [
        (
            "Upstream_intergenic",
            lambda: ref_gene.getIntergenic(direction="up", size=upstream_size),
        ),
        ("Five_prime_UTR", lambda: ref_gene.getUTRs(utr=5)),
        ("Coding_exon", ref_gene.getCDSExons),
        ("First_intron", lambda: ref_gene.getIntrons(itype="first")),
        ("Internal_intron", lambda: ref_gene.getIntrons(itype="internal")),
        ("Last_intron", lambda: ref_gene.getIntrons(itype="last")),
        ("Three_prime_UTR", lambda: ref_gene.getUTRs(utr=3)),
        (
            "Downstream_intergenic",
            lambda: ref_gene.getIntergenic(direction="down", size=downstream_size),
        ),
    ]

    profiles = []
    for region_name, builder in region_builders:
        printlog(f"Processing {REGION_LABELS[region_name]} regions ...")
        profiles.append(
            compute_region_profile(region_name, builder(), cpg_ranges)
        )

    nonempty = [profile for profile in profiles if not profile.empty]
    if not nonempty:
        raise GeneProfileError(
            "No methylation profile could be generated for any genomic region."
        )

    return pd.concat(nonempty, ignore_index=True)


def plot_profile(profile: pd.DataFrame, output_path: Path, dpi: int, width: float, height: float) -> None:
    """Plot concatenated profiles with boundaries between genomic regions."""
    ordered_groups = [name for name in REGION_LABELS if name in set(profile["Group"])]

    x_values: list[int] = []
    y_values: list[float] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    boundaries: list[int] = []
    cursor = 0

    for group_name in ordered_groups:
        group = profile.loc[profile["Group"] == group_name].sort_values(
            "Relative_position(5'->3')"
        )
        values = group["Average_beta"].to_numpy(dtype=float)
        if values.size == 0:
            continue

        start = cursor
        x_values.extend(range(cursor, cursor + values.size))
        y_values.extend(values.tolist())
        cursor += values.size

        tick_positions.append((start + cursor - 1) / 2.0)
        tick_labels.append(REGION_LABELS[group_name])
        boundaries.append(cursor - 0.5)

    if not y_values:
        raise GeneProfileError("No finite values are available for plotting.")

    figure, axis = plt.subplots(figsize=(width, height))
    axis.plot(x_values, y_values, linewidth=1.5)
    axis.axhline(0.5, linestyle="--", linewidth=0.8)

    for boundary in boundaries[:-1]:
        axis.axvline(boundary, linestyle="--", linewidth=0.8)

    axis.set_xlim(min(x_values), max(x_values))
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("Average methylation")
    axis.set_xlabel("Transcript-oriented genomic regions (5' to 3')")
    axis.set_xticks(tick_positions)
    axis.set_xticklabels(tick_labels, rotation=30, ha="right")
    axis.set_title("Gene-centered DNA methylation profile")
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
    table_path = Path(f"{args.out_prefix}.gene_centered_profile.tsv")

    try:
        profile = generate_profiles(
            input_file=args.input_file,
            gene_file=args.refgene,
            upstream_size=args.upstream,
            downstream_size=args.downstream,
        )

        printlog(f'Writing profile table: "{table_path}"')
        profile.to_csv(table_path, sep="\t", index=False, float_format="%.6f")

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)
            for extension in extensions:
                plot_path = Path(f"{args.out_prefix}.gene_centered_profile.{extension}")
                printlog(f'Writing profile plot: "{plot_path}"')
                plot_profile(
                    profile=profile,
                    output_path=plot_path,
                    dpi=args.dpi,
                    width=args.width,
                    height=args.height,
                )
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

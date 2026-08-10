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

"""Summarize CpG distribution across gene-centered genomic annotations.

The program partitions the genome into five prioritized annotation classes:

0. Coding exons
1. UTR exons
2. Introns
3. Upstream regions of transcription start sites (TSS)
4. Downstream regions of transcription end sites (TES)

Higher-priority annotations override lower-priority annotations. For example,
bases annotated as both exon and intron are assigned to the exon category.

The output reports, for each category:

* number of non-overlapping genomic intervals
* total annotated bases
* raw CpG count
* CpG density per kilobase

The reference gene model must be BED12.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd

from cpgmodule import BED
from cpgmodule._version import __version__
from cpgmodule.utils import count_over_range, printlog, read_CpG_bed


class GeneCenteredDistributionError(RuntimeError):
    """Raised when gene-centered CpG distribution analysis fails."""


PRIORITY = (
    ("0", "Coding exons"),
    ("1", "UTR exons"),
    ("2", "Introns"),
    ("3", "Upstream of TSS"),
    ("4", "Downstream of TES"),
)


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools CpG_distrb_gene_centered",
        description=(
            "Calculate CpG distribution across prioritized gene-centered "
            "genomic annotations."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input CpG BED3+ file. Only chromosome/start/end are required. "
            "Compressed input is supported by the CpGtools reader."
        ),
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="gene_file",
        required=True,
        type=Path,
        help="Reference gene model in BED12 format.",
    )
    parser.add_argument(
        "-u",
        "--upstream",
        dest="upstream_size",
        type=int,
        default=2000,
        help="Size of the upstream region relative to the TSS, in bp.",
    )
    parser.add_argument(
        "-d",
        "--downstream",
        dest="downstream_size",
        type=int,
        default=2000,
        help="Size of the downstream region relative to the TES, in bp.",
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
        default=8.0,
        help="Plot width in inches.",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=6.0,
        help="Plot height in inches.",
    )
    parser.add_argument(
        "--no_plot",
        action="store_true",
        help="Write the summary table without generating a plot.",
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

    if args.dpi < 1:
        parser.error("--dpi must be at least 1")

    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def safe_union(intervals):
    """Merge BED3 intervals and return an empty list when input is empty."""
    if not intervals:
        return []
    return BED.unionBed3(intervals)


def safe_subtract(intervals, higher_priority):
    """Subtract higher-priority BED3 intervals when both collections are nonempty."""
    if not intervals:
        return []
    if not higher_priority:
        return intervals
    return BED.subtractBed3(intervals, higher_priority)


def summarize_category(name: str, intervals, cpg_ranges) -> dict[str, object]:
    """Calculate interval count, total size, CpG count, and CpG density."""
    if not intervals:
        return {
            "Name": name,
            "Number_of_regions": 0,
            "Size_of_regions_bp": 0,
            "CpG_raw_count": 0,
            "CpG_count_per_KB": 0.0,
        }

    size, count = count_over_range(intervals, cpg_ranges)

    density = 0.0
    if size > 0:
        density = count * 1000.0 / size

    return {
        "Name": name,
        "Number_of_regions": len(intervals),
        "Size_of_regions_bp": int(size),
        "CpG_raw_count": int(count),
        "CpG_count_per_KB": float(density),
    }


def build_prioritized_regions(
    ref_gene: BED.ParseBED,
    upstream_size: int,
    downstream_size: int,
):
    """Construct non-overlapping annotation classes in priority order."""
    printlog("Extracting coding exons ...")
    coding = ref_gene.getCDSExons(
        uniquify=True,
        stranded=False,
    )
    coding = safe_union(coding)

    printlog("Extracting UTR exons ...")
    utr = ref_gene.getUTRs(
        utr=35,
        uniquify=True,
        stranded=False,
    )
    utr = safe_union(utr)
    utr = safe_subtract(utr, coding)

    printlog("Extracting introns ...")
    introns = ref_gene.getIntrons(
        itype="all",
        uniquify=True,
        stranded=False,
    )
    introns = safe_union(introns)
    introns = safe_subtract(introns, coding)
    introns = safe_subtract(introns, utr)

    printlog("Extracting upstream regions ...")
    upstream = ref_gene.getIntergenic(
        direction="up",
        size=upstream_size,
        uniquify=True,
        stranded=False,
    )
    upstream = safe_union(upstream)
    upstream = safe_subtract(upstream, coding)
    upstream = safe_subtract(upstream, utr)
    upstream = safe_subtract(upstream, introns)

    printlog("Extracting downstream regions ...")
    downstream = ref_gene.getIntergenic(
        direction="down",
        size=downstream_size,
        uniquify=True,
        stranded=False,
    )
    downstream = safe_union(downstream)
    downstream = safe_subtract(downstream, coding)
    downstream = safe_subtract(downstream, utr)
    downstream = safe_subtract(downstream, introns)
    downstream = safe_subtract(downstream, upstream)

    return {
        "Coding exons": coding,
        "UTR exons": utr,
        "Introns": introns,
        "Upstream of TSS": upstream,
        "Downstream of TES": downstream,
    }


def generate_summary(
    input_file: Path,
    gene_file: Path,
    upstream_size: int,
    downstream_size: int,
) -> pd.DataFrame:
    """Generate the prioritized gene-centered CpG distribution table."""
    printlog(f'Reading CpG file: "{input_file}"')
    cpg_ranges = read_CpG_bed(str(input_file))

    printlog(f'Reading reference gene model: "{gene_file}"')
    ref_gene = BED.ParseBED(str(gene_file))

    regions = build_prioritized_regions(
        ref_gene=ref_gene,
        upstream_size=upstream_size,
        downstream_size=downstream_size,
    )

    rows = []

    for priority_order, name in PRIORITY:
        printlog(f'Counting CpGs in "{name}" ...')
        summary = summarize_category(
            name=name,
            intervals=regions[name],
            cpg_ranges=cpg_ranges,
        )
        summary["Priority_order"] = priority_order
        rows.append(summary)

    result = pd.DataFrame(rows)
    return result[
        [
            "Priority_order",
            "Name",
            "Number_of_regions",
            "Size_of_regions_bp",
            "CpG_raw_count",
            "CpG_count_per_KB",
        ]
    ]


def plot_summary(
    summary: pd.DataFrame,
    output_path: Path,
    *,
    dpi: int,
    width: float,
    height: float,
) -> None:
    """Create a bar plot of CpG density by genomic annotation class."""
    figure, axis = plt.subplots(figsize=(width, height))

    x = range(len(summary))
    axis.bar(
        x,
        summary["CpG_count_per_KB"].to_numpy(dtype=float),
        width=0.75,
    )

    axis.set_xticks(list(x))
    axis.set_xticklabels(
        summary["Name"].tolist(),
        rotation=30,
        ha="right",
    )
    axis.set_ylabel("CpG count per kb")
    axis.set_xlabel("Genomic annotation")
    axis.set_title("CpG density across gene-centered genomic regions")

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
    table_path = Path(f"{args.out_prefix}.gene_centered_distribution.tsv")

    try:
        summary = generate_summary(
            input_file=args.input_file,
            gene_file=args.gene_file,
            upstream_size=args.upstream_size,
            downstream_size=args.downstream_size,
        )

        printlog(f'Writing CpG distribution table: "{table_path}"')
        summary.to_csv(
            table_path,
            sep="\t",
            index=False,
            float_format="%.6f",
        )

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)

            for extension in extensions:
                plot_path = Path(
                    f"{args.out_prefix}.gene_centered_distribution.{extension}"
                )
                printlog(f'Writing CpG distribution plot: "{plot_path}"')
                plot_summary(
                    summary=summary,
                    output_path=plot_path,
                    dpi=args.dpi,
                    width=args.width,
                    height=args.height,
                )

    except GeneCenteredDistributionError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

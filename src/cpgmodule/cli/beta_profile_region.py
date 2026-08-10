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

"""Generate transcript-oriented DNA methylation profiles around genomic regions.

The program calculates average methylation across three normalized region classes:

* upstream flanking regions
* user-supplied genomic regions
* downstream flanking regions

The methylation input must be BED6 or BED6+ with beta values in column 5.
The target-region input must be BED3 or BED6+. When strand is absent, regions are
assumed to be on the ``+`` strand. All profiles are oriented 5' to 3' relative to
the target-region strand.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd

from cpgmodule._version import __version__
from cpgmodule.utils import (
    coverage_over_range,
    printlog,
    read_CpG_bed,
    read_region_bed,
)


class RegionProfileError(RuntimeError):
    """Raised when profile generation cannot be completed."""


PROFILE_LABELS = {
    "Upstream_region": "Upstream",
    "User_region": "Target region",
    "Downstream_region": "Downstream",
}


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_profile_region",
        description=(
            "Generate a strand-aware DNA methylation profile across user-defined "
            "regions and their upstream/downstream flanks."
        ),
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
        "--region",
        required=True,
        type=Path,
        dest="region_file",
        help=(
            "BED3+ file containing target regions. If column 6 is absent, the "
            "region is treated as '+' strand."
        ),
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=2000,
        help="Length of the upstream flank in base pairs.",
    )
    parser.add_argument(
        "-d",
        "--downstream",
        type=int,
        default=2000,
        help="Length of the downstream flank in base pairs.",
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        "--output",
        required=True,
        type=Path,
        dest="out_prefix",
        help="Output filename prefix, optionally including a directory.",
    )
    parser.add_argument(
        "--format",
        choices=("png", "pdf", "both"),
        default="pdf",
        help="Plot output format.",
    )
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution.")
    parser.add_argument("--width", type=float, default=8.0, help="Plot width in inches.")
    parser.add_argument("--height", type=float, default=5.0, help="Plot height in inches.")
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
    if not args.region_file.is_file():
        parser.error(f"Region file does not exist: {args.region_file}")
    if args.upstream < 0:
        parser.error("--upstream must be non-negative")
    if args.downstream < 0:
        parser.error("--downstream must be non-negative")
    if args.dpi < 1:
        parser.error("--dpi must be at least 1")
    if args.width <= 0 or args.height <= 0:
        parser.error("--width and --height must be greater than zero")


def load_regions(path: Path) -> list[tuple[str, int, int, str]]:
    """Read, validate, and deduplicate target regions."""
    regions: list[tuple[str, int, int, str]] = []
    for chrom, start, end, strand in read_region_bed(str(path)):
        start = int(start)
        end = int(end)
        strand = strand if strand in {"+", "-"} else "+"

        if start < 0:
            raise RegionProfileError(
                f"Negative start coordinate in region file: {chrom}:{start}-{end}"
            )
        if end <= start:
            raise RegionProfileError(
                f"Region end must be greater than start: {chrom}:{start}-{end}"
            )

        regions.append((str(chrom), start, end, strand))

    if not regions:
        raise RegionProfileError(f"No valid regions were found in {path}")

    # Preserve file order while removing exact duplicate intervals.
    return list(dict.fromkeys(regions))


def build_flanks(
    regions: Iterable[tuple[str, int, int, str]],
    size: int,
    direction: str,
) -> list[tuple[str, int, int, str]]:
    """Construct strand-aware upstream or downstream flanking intervals."""
    if direction not in {"upstream", "downstream"}:
        raise ValueError("direction must be 'upstream' or 'downstream'")
    if size <= 0:
        return []

    flanks: list[tuple[str, int, int, str]] = []
    for chrom, start, end, strand in regions:
        if direction == "upstream":
            if strand == "+":
                flank_start = max(0, start - size)
                flank_end = start
            else:
                flank_start = end
                flank_end = end + size
        else:
            if strand == "+":
                flank_start = end
                flank_end = end + size
            else:
                flank_start = max(0, start - size)
                flank_end = start

        if flank_end > flank_start:
            flanks.append((chrom, flank_start, flank_end, strand))

    return list(dict.fromkeys(flanks))


def compute_profile(
    group_name: str,
    regions: list[tuple[str, int, int, str]],
    cpg_ranges,
) -> pd.DataFrame:
    """Calculate a normalized methylation profile for one region group."""
    if not regions:
        return pd.DataFrame(
            columns=["Group", "Relative_position(5'->3')", "Average_beta"]
        )

    values = coverage_over_range(regions, cpg_ranges)
    if not values:
        printlog(f'Warning: no methylation values were found for "{group_name}".')
        return pd.DataFrame(
            columns=["Group", "Relative_position(5'->3')", "Average_beta"]
        )

    rows = [
        {
            "Group": group_name,
            "Relative_position(5'->3')": int(position),
            "Average_beta": float(values[position]),
        }
        for position in sorted(values)
    ]
    return pd.DataFrame(rows)


def generate_profile(
    input_file: Path,
    region_file: Path,
    upstream_size: int,
    downstream_size: int,
) -> pd.DataFrame:
    """Generate upstream, target-region, and downstream methylation profiles."""
    printlog(f'Reading CpG file: "{input_file}"')
    cpg_ranges = read_CpG_bed(str(input_file))

    printlog(f'Reading region file: "{region_file}"')
    target_regions = load_regions(region_file)
    printlog(f"Loaded {len(target_regions)} unique target regions")

    upstream_regions = build_flanks(target_regions, upstream_size, "upstream")
    downstream_regions = build_flanks(target_regions, downstream_size, "downstream")

    region_groups = [
        ("Upstream_region", upstream_regions),
        ("User_region", target_regions),
        ("Downstream_region", downstream_regions),
    ]

    profiles: list[pd.DataFrame] = []
    for group_name, regions in region_groups:
        printlog(f"Calculating {PROFILE_LABELS[group_name]} profile ...")
        profile = compute_profile(group_name, regions, cpg_ranges)
        if not profile.empty:
            profiles.append(profile)

    if not profiles:
        raise RegionProfileError(
            "No methylation profile could be generated for the target or flanking regions."
        )

    return pd.concat(profiles, ignore_index=True)


def plot_profile(
    profile: pd.DataFrame,
    output_path: Path,
    dpi: int,
    width: float,
    height: float,
) -> None:
    """Plot concatenated upstream, target, and downstream profiles."""
    ordered_groups = [
        group for group in PROFILE_LABELS if group in set(profile["Group"])
    ]

    x_values: list[int] = []
    y_values: list[float] = []
    tick_positions: list[float] = []
    tick_labels: list[str] = []
    boundaries: list[float] = []
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
        tick_labels.append(PROFILE_LABELS[group_name])
        boundaries.append(cursor - 0.5)

    if not y_values:
        raise RegionProfileError("No finite values are available for plotting.")

    figure, axis = plt.subplots(figsize=(width, height))
    axis.plot(x_values, y_values, linewidth=1.5)
    axis.axhline(0.5, linestyle="--", linewidth=0.8)

    for boundary in boundaries[:-1]:
        axis.axvline(boundary, linestyle="--", linewidth=0.8)

    axis.set_xlim(min(x_values), max(x_values))
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel("Average methylation")
    axis.set_xlabel("Transcript-oriented regions (5' to 3')")
    axis.set_xticks(tick_positions)
    axis.set_xticklabels(tick_labels, rotation=25, ha="right")
    axis.set_title("Region-centered DNA methylation profile")
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
    table_path = Path(f"{args.out_prefix}.region_profile.tsv")

    try:
        profile = generate_profile(
            input_file=args.input_file,
            region_file=args.region_file,
            upstream_size=args.upstream,
            downstream_size=args.downstream,
        )

        printlog(f'Writing profile table: "{table_path}"')
        profile.to_csv(table_path, sep="\t", index=False, float_format="%.6f")

        if not args.no_plot:
            extensions = ("png", "pdf") if args.format == "both" else (args.format,)
            for extension in extensions:
                plot_path = Path(f"{args.out_prefix}.region_profile.{extension}")
                printlog(f'Writing profile plot: "{plot_path}"')
                plot_profile(
                    profile=profile,
                    output_path=plot_path,
                    dpi=args.dpi,
                    width=args.width,
                    height=args.height,
                )
    except (RegionProfileError, OSError, ValueError) as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

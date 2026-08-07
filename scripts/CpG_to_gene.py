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
Annotate CpGs by assigning them to putative target genes using the
"basal plus extension" rules used by GREAT.

Basal regulatory domain
~~~~~~~~~~~~~~~~~~~~~~~
A user-defined genomic region around the transcription start site (TSS).
By default, the basal domain extends 5 kb upstream and 1 kb downstream of
the TSS. Nearby genes are ignored when defining a gene's basal domain, so
basal domains from different genes may overlap.

Extended regulatory domain
~~~~~~~~~~~~~~~~~~~~~~~~~~
A gene's regulatory domain is extended in both directions to the nearest
gene's basal regulatory domain, but by no more than the maximum extension
distance (default: 1,000 kb) in either direction.

Notes
-----
1. CpG-to-gene assignment depends strongly on the gene annotation used.
   A conservative gene model, such as curated protein-coding RefSeq genes,
   is recommended.
2. Multiple isoforms of the same gene should ideally be merged into one
   representative gene model before running this command.
"""

import argparse
import sys
from pathlib import Path

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.region2gene import getBasalDomains, geteExtendedDomains
from cpgmodule.utils import printlog


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
            "BED3+ file specifying CpG positions. The file may be plain text "
            "or compressed (.gz, .bz2)."
        ),
    )

    parser.add_argument(
        "-r",
        "--refgene",
        dest="gene_file",
        required=True,
        help=(
            "Reference gene model in BED12 format. One gene per transcript "
            "record is recommended; for genes with multiple isoforms, use a "
            "collapsed or canonical transcript representation."
        ),
    )

    parser.add_argument(
        "-u",
        "--basal-up",
        dest="basal_up_size",
        type=int,
        default=5000,
        help=(
            "Upstream extension from the TSS used to define the basal "
            "regulatory domain [default: %(default)s bp]."
        ),
    )

    parser.add_argument(
        "-d",
        "--basal-down",
        dest="basal_down_size",
        type=int,
        default=1000,
        help=(
            "Downstream extension from the TSS used to define the basal "
            "regulatory domain [default: %(default)s bp]."
        ),
    )

    parser.add_argument(
        "-e",
        "--extension",
        dest="extension_size",
        type=int,
        default=1_000_000,
        help=(
            "Maximum extension from the TSS used to define the extended "
            "regulatory domain [default: %(default)s bp]."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help=(
            "Output prefix. Results are written to "
            "<prefix>.associated_genes.txt."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(args, parser=None):
    """Validate parsed command-line arguments."""
    if args.basal_up_size < 0:
        message = "--basal-up must be zero or greater"
        if parser is not None:
            parser.error(message)
        raise ValueError(message)

    if args.basal_down_size < 0:
        message = "--basal-down must be zero or greater"
        if parser is not None:
            parser.error(message)
        raise ValueError(message)

    if args.extension_size < 0:
        message = "--extension must be zero or greater"
        if parser is not None:
            parser.error(message)
        raise ValueError(message)


def parse_bed_record(line, line_number):
    """
    Parse a BED3+ record.

    Returns
    -------
    tuple or None
        (chrom, start, end) for a valid record, otherwise None.
    """
    fields = line.split()

    if len(fields) < 3:
        print(
            f"Invalid BED line {line_number}: expected at least 3 columns: "
            f"{line}",
            file=sys.stderr,
        )
        return None

    try:
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
    except ValueError:
        print(
            f"Invalid BED coordinates on line {line_number}: {line}",
            file=sys.stderr,
        )
        return None

    if start < 0:
        print(
            f"Invalid BED start on line {line_number}: {start}",
            file=sys.stderr,
        )
        return None

    if end < start:
        print(
            f"BED start cannot be greater than end on line {line_number}: "
            f"{line}",
            file=sys.stderr,
        )
        return None

    return chrom, start, end


def overlapping_gene_names(domain_map, chrom, start, end):
    """
    Return gene names whose regulatory domains overlap a BED interval.

    The historical CpGtools sentinel "//" is returned when no gene overlaps.
    """
    if chrom not in domain_map:
        return {"//"}

    overlaps = domain_map[chrom].find(start, end)
    if not overlaps:
        return {"//"}

    return {overlap.value for overlap in overlaps}


def annotate_cpgs(
    input_file,
    gene_file,
    output_prefix,
    basal_up_size=5000,
    basal_down_size=1000,
    extension_size=1_000_000,
):
    """
    Annotate CpGs with genes whose basal or extended domains overlap them.

    Returns
    -------
    str
        Path to the generated output file.
    """
    output_path = Path(f"{output_prefix}.associated_genes.txt")

    if output_path.parent != Path("."):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    printlog(
        f'Calculate basal regulatory domain from: "{gene_file}" ...'
    )
    basal_domains = getBasalDomains(
        bedfile=gene_file,
        up=basal_up_size,
        down=basal_down_size,
        printit=False,
    )

    printlog(
        f'Calculate extended regulatory domain from: "{gene_file}" ...'
    )
    extended_domains = geteExtendedDomains(
        basal_ranges=basal_domains,
        bedfile=gene_file,
        up=basal_up_size,
        down=basal_down_size,
        ext=extension_size,
        printit=False,
    )

    printlog("Assigning CpG to gene ...")

    with open(output_path, "w") as fout:
        print(
            "#The last column contains genes whose extended regulatory "
            "domain are overlapped with the CpG",
            file=fout,
        )
        print(
            "#The 2nd last column contains genes whose basal regulatory "
            "domain are overlapped with the CpG",
            file=fout,
        )
        print('#"//" indicates no genes are found', file=fout)

        for line_number, raw_line in enumerate(
            ireader.reader(input_file), start=1
        ):
            line = raw_line.rstrip("\r\n")

            if not line:
                continue

            if line.startswith("#"):
                print(line, file=fout)
                continue

            if line.startswith("track") or line.startswith("browser"):
                continue

            parsed = parse_bed_record(line, line_number)
            if parsed is None:
                continue

            chrom, start, end = parsed

            basal_genes = overlapping_gene_names(
                basal_domains,
                chrom,
                start,
                end,
            )

            extended_genes = overlapping_gene_names(
                extended_domains,
                chrom,
                start,
                end,
            )

            # Preserve historical behavior: genes already reported in the
            # basal column are removed from the extended-only column.
            extended_only = extended_genes - basal_genes

            if not extended_only:
                extended_only = {"//"}

            print(
                line
                + "\t"
                + ";".join(sorted(basal_genes))
                + "\t"
                + ";".join(sorted(extended_only)),
                file=fout,
            )

    return str(output_path)


def main(argv=None):
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str] or None
        Optional argument list. When None, argparse reads sys.argv.
        Passing a list allows this function to be called cleanly from a
        CpGtools dispatcher or from tests.

    Returns
    -------
    int
        Process-style return code: 0 on success, nonzero on failure.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        validate_args(args, parser=parser)

        annotate_cpgs(
            input_file=args.input_file,
            gene_file=args.gene_file,
            output_prefix=args.out_file,
            basal_up_size=args.basal_up_size,
            basal_down_size=args.basal_down_size,
            extension_size=args.extension_size,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

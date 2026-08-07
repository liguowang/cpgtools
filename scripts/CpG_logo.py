#!/usr/bin/env python3

"""
Description
-----------
Generate a DNA motif logo and motif matrices for a set of CpG positions.

The input BED file should contain at least six columns:

    chrom  chromStart  chromEnd  name  score  strand

The third BED column (chromEnd) is treated as the genomic position of the
methylated cytosine, preserving the historical behavior of CpGtools.

For each valid CpG record, this script:
  1. extracts a strand-aware sequence window from the reference FASTA;
  2. writes the extracted sequences to FASTA;
  3. generates WebLogo PDF and PNG files when `weblogo` is available;
  4. writes PFM, PPM, PWM, JASPAR, and MEME motif files.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pysam

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.imotif import PSSM
from cpgmodule.utils import printlog, revcomp


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
            "BED file specifying CpG positions. At least six columns are "
            "expected: chrom, chromStart, chromEnd, name, score, strand. "
            "Compressed .gz and .bz2 files are supported by cpgmodule.ireader."
        ),
    )
    parser.add_argument(
        "-r",
        "--refgenome",
        dest="genome_file",
        required=True,
        help=(
            "Reference genome in FASTA format. If the corresponding .fai "
            "index is missing, it will be created with pysam.faidx()."
        ),
    )
    parser.add_argument(
        "-e",
        "--extend",
        dest="extend_size",
        type=int,
        default=5,
        help=(
            "Number of bases to extend upstream and downstream "
            "[default: %(default)s bp]."
        ),
    )
    parser.add_argument(
        "-n",
        "--name",
        dest="motif_name",
        default="motif",
        help="Motif name [default: %(default)s].",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help="Prefix for output files.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def ensure_fasta_index(genome_file):
    """Create a FASTA index when it does not already exist."""
    genome_path = Path(genome_file)

    if not genome_path.is_file():
        raise FileNotFoundError(
            f"Reference genome does not exist: {genome_file}"
        )

    index_path = Path(f"{genome_file}.fai")

    if not index_path.exists():
        printlog(f"Creating index for {genome_file}")
        pysam.faidx(str(genome_path))


def parse_bed_record(line, line_number):
    """
    Parse one BED record.

    Returns
    -------
    tuple or None
        (chrom, position, strand) for a valid record, otherwise None.
    """
    fields = line.split()

    if len(fields) < 6:
        print(
            f"BED must have at least 6 columns. "
            f"Skip line {line_number}: {line}",
            file=sys.stderr,
        )
        return None

    chrom = fields[0]

    try:
        position = int(fields[2])
    except ValueError:
        print(
            f"Invalid BED coordinate in column 3. "
            f"Skip line {line_number}: {line}",
            file=sys.stderr,
        )
        return None

    strand = fields[5]
    if strand not in {"+", "-"}:
        print(
            f"Invalid strand '{strand}'. Expected '+' or '-'. "
            f"Skip line {line_number}: {line}",
            file=sys.stderr,
        )
        return None

    return chrom, position, strand


def extract_sequences(input_file, genome_file, extend_size, fasta_output):
    """
    Extract strand-aware sequence windows around CpG positions.

    The coordinate calculation intentionally preserves the original CpGtools
    behavior:

        start = chromEnd - extend_size - 1
        end   = chromEnd + extend_size

    Therefore, each sequence has length 2 * extend_size + 1.

    Returns
    -------
    int
        Number of sequences written.
    """
    sequence_count = 0
    ambiguous_count = 0

    printlog(f"Reading {input_file} ...")

    with pysam.FastaFile(genome_file) as reference, open(
        fasta_output, "w"
    ) as fout:
        reference_names = set(reference.references)

        for line_number, raw_line in enumerate(
            ireader.reader(input_file), start=1
        ):
            line = raw_line.strip()

            if not line:
                continue

            if (
                line.startswith("#")
                or line.startswith("track")
                or line.startswith("browser")
            ):
                continue

            parsed = parse_bed_record(line, line_number)
            if parsed is None:
                continue

            chrom, position, strand = parsed

            if chrom not in reference_names:
                print(
                    f"Chromosome '{chrom}' is not present in the reference "
                    f"genome. Skip line {line_number}.",
                    file=sys.stderr,
                )
                continue

            start = position - extend_size - 1
            end = position + extend_size

            if start < 0 or start >= end:
                print(
                    f"Invalid interval {chrom}:{start}-{end}. "
                    f"Skip line {line_number}.",
                    file=sys.stderr,
                )
                continue

            chrom_length = reference.get_reference_length(chrom)
            if end > chrom_length:
                print(
                    f"Interval {chrom}:{start}-{end} exceeds chromosome "
                    f"length ({chrom_length}). Skip line {line_number}.",
                    file=sys.stderr,
                )
                continue

            try:
                sequence = reference.fetch(chrom, start, end).upper()
            except (KeyError, ValueError) as exc:
                print(
                    f"Cannot fetch {chrom}:{start}-{end} on line "
                    f"{line_number}: {exc}",
                    file=sys.stderr,
                )
                continue

            if strand == "-":
                sequence = revcomp(sequence)

            # cpgmodule.imotif.PSSM accepts DNA sequences containing only
            # A/C/G/T. Skip windows containing N or other ambiguous symbols
            # here so downstream matrix generation stays clean and consistent.
            if set(sequence) - {"A", "C", "G", "T"}:
                ambiguous_count += 1
                continue

            fasta_name = f">{chrom}_{start}_{end}_{strand}"
            fout.write(f"{fasta_name}\n{sequence}\n")
            sequence_count += 1

    if ambiguous_count:
        printlog(
            f"Skipped {ambiguous_count} sequence(s) containing ambiguous "
            "reference bases (for example, N)"
        )

    return sequence_count


def generate_weblogo(fasta_file, output_prefix, motif_name):
    """
    Generate PDF and PNG sequence logos.

    WebLogo requires Ghostscript for PDF and PNG output. Detect that dependency
    before invoking WebLogo so a missing `gs` executable produces a concise
    warning rather than a Python traceback from WebLogo.

    Returns
    -------
    bool
        True when both logos were generated successfully, otherwise False.
    """
    weblogo = shutil.which("weblogo")

    if weblogo is None:
        print(
            "Warning: 'weblogo' was not found in PATH. "
            "Skipping logo generation.",
            file=sys.stderr,
        )
        return False

    ghostscript = shutil.which("gs")
    if ghostscript is None:
        print(
            "Warning: Ghostscript ('gs') was not found in PATH. "
            "WebLogo requires Ghostscript for PDF/PNG output, so logo "
            "generation is being skipped. Motif matrices will still be "
            "generated.\n"
            "Install Ghostscript and rerun CpG_logo.py to create the logos.",
            file=sys.stderr,
        )
        return False

    output_formats = (
        ("PDF", f"{output_prefix}.logo.pdf"),
        ("PNG", f"{output_prefix}.logo.png"),
    )

    all_ok = True

    for output_format, output_file in output_formats:
        command = [
            weblogo,
            "--format",
            output_format,
            "-D",
            "fasta",
            "-c",
            "classic",
            "-s",
            "large",
            "-f",
            str(fasta_file),
            "-o",
            output_file,
            "-t",
            motif_name,
        ]

        try:
            result = subprocess.run(
                command,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except subprocess.CalledProcessError as exc:
            all_ok = False
            detail = (exc.stderr or "").strip()
            print(
                f"Warning: WebLogo failed while generating {output_format} "
                f"output: {output_file}",
                file=sys.stderr,
            )
            if detail:
                # Show WebLogo's useful diagnostic without a full subprocess
                # traceback from this wrapper.
                print(detail, file=sys.stderr)

    if all_ok:
        printlog(
            f'Motif logo saved to "{output_prefix}.logo.pdf" and '
            f'"{output_prefix}.logo.png"'
        )

    return all_ok

def write_motif_matrices(fasta_file, output_prefix, motif_name):
    """Write PFM, PPM, PWM, JASPAR, and MEME motif files."""
    motif = PSSM(sites=str(fasta_file), name=motif_name)

    outputs = (
        (".pfm", "position frequency matrix (PFM)", motif.toPFM),
        (".ppm", "position probability matrix (PPM)", motif.toPPM),
        (".pwm", "position weight matrix (PWM)", motif.toPWM),
        (".jaspar", "JASPAR format matrix", motif.toJaspar),
        (".meme", "MEME format matrix", motif.toMEME),
    )

    for suffix, description, writer in outputs:
        output_file = f"{output_prefix}{suffix}"
        printlog(f'Write {description} to "{output_file}"')

        with open(output_file, "w") as fout:
            writer(FOUT=fout)


def run(
    input_file,
    genome_file,
    output_prefix,
    extend_size=5,
    motif_name="motif",
):
    """
    Run CpG sequence extraction, logo generation, and motif export.

    This function contains the reusable application logic independently of
    command-line parsing, making it easier to test from within CpGtools.
    """
    if extend_size < 0:
        raise ValueError("extend_size must be zero or greater")

    ensure_fasta_index(genome_file)

    output_prefix = Path(output_prefix)

    if output_prefix.parent != Path("."):
        output_prefix.parent.mkdir(parents=True, exist_ok=True)

    fasta_output = Path(f"{output_prefix}.fa")

    sequence_count = extract_sequences(
        input_file=input_file,
        genome_file=genome_file,
        extend_size=extend_size,
        fasta_output=fasta_output,
    )

    if sequence_count == 0:
        raise RuntimeError(
            "No valid sequences were extracted from the input BED file."
        )

    printlog(
        f'Extracted {sequence_count} sequence(s) to "{fasta_output}"'
    )

    printlog("Generate motif logo ...")
    generate_weblogo(
        fasta_file=fasta_output,
        output_prefix=str(output_prefix),
        motif_name=motif_name,
    )

    write_motif_matrices(
        fasta_file=fasta_output,
        output_prefix=str(output_prefix),
        motif_name=motif_name,
    )

    return sequence_count


def main(argv=None):
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        run(
            input_file=args.input_file,
            genome_file=args.genome_file,
            output_prefix=args.out_file,
            extend_size=args.extend_size,
            motif_name=args.motif_name,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

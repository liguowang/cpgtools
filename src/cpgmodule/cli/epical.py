#!/usr/bin/env python3
"""Command-line interface for CpGtools epigenetic clocks."""

import argparse
import sys

from cpgmodule._version import __version__
from dmc import helpdoc
from dmc.clock_info import clockinfo

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


# Map CLI clock names to the model metadata used for help text.
CLOCK_MODELS = {
    "Horvath13": "Horvath13.pkl",
    "Horvath13_shrunk": "Horvath13_shrunk.pkl",
    "Horvath18": "Horvath18.pkl",
    "Levine": "Levine.pkl",
    "Hannum": "Hannum.pkl",
    "Zhang_EN": "Zhang_EN.pkl",
    "Zhang_BLUP": "Zhang_BLUP.pkl",
    "AltumAge": "AltumAge.pkl",
    "Lu_DNAmTL": "Lu_DNAmTL.pkl",
    "Weidner": "Weidner.pkl",
    "Lin": "Lin.pkl",
    "ENCen100": "ENCen100.pkl",
    "ENCen40": "ENCen40.pkl",
    "DunedinPACE": "DunedinPACE.pkl",
    "Ped_Wu": "Ped_Wu.pkl",
    "PedBE": "PedBE.pkl",
    "GA_Bohlin": "GA_Bohlin.pkl",
    "GA_Haftorn": "GA_Haftorn.pkl",
    "GA_Knight": "GA_Knight.pkl",
    "GA_Mayne": "GA_Mayne.pkl",
    "GA_Lee_CPC": "GA_Lee_CPC.pkl",
    "GA_Lee_RPC": "GA_Lee_RPC.pkl",
    "GA_Lee_rRPC": "GA_Lee_rRPC.pkl",
    "Cortical": "Cortical.pkl",
    "MEAT": "MEAT.pkl",
    "WLMT": "WLMT_mm10.pkl",
    "YOMT": "YOMT_mm10.pkl",
    "mmLiver": "mmLiver_mm10.pkl",
    "mmBlood": "mmBlood_mm10.pkl",
    "mammClock1": "mammClock1.pkl",
    "mammClock2": "mammClock2.pkl",
    "mammClock3": "mammClock3.pkl",
}

STANDARD_CLOCKS = (
    "Horvath13",
    "Horvath13_shrunk",
    "Horvath18",
    "Levine",
    "Hannum",
    "Zhang_EN",
    "Zhang_BLUP",
    "AltumAge",
    "Lu_DNAmTL",
    "Ped_Wu",
    "PedBE",
    "GA_Bohlin",
    "GA_Haftorn",
    "GA_Knight",
    "GA_Mayne",
    "GA_Lee_CPC",
    "GA_Lee_RPC",
    "GA_Lee_rRPC",
    "Cortical",
    "MEAT",
    "Weidner",
    "Lin",
    "ENCen100",
    "ENCen40",
    "mammClock1",
)

MOUSE_CLOCKS = ("WLMT", "YOMT", "mmLiver", "mmBlood")
MAMMALIAN_CLOCKS = ("mammClock2", "mammClock3")

HORVATH_CLOCKS = {
    "Horvath13",
    "Horvath13_shrunk",
    "Horvath18",
    "MEAT",
    "PedBE",
    "Cortical",
}
LEVINE_HANNUM_CLOCKS = {"Levine", "Hannum", "Lu_DNAmTL"}
ZHANG_CLOCKS = {"Zhang_BLUP", "Zhang_EN"}
GESTATIONAL_CLOCKS = {
    "GA_Knight",
    "GA_Mayne",
    "GA_Bohlin",
    "GA_Haftorn",
    "GA_Lee_CPC",
    "GA_Lee_RPC",
    "GA_Lee_rRPC",
}
GENERAL_CLOCKS = {"Weidner", "Lin", "ENCen100", "ENCen40"}


def _clock_help(name):
    """Return the existing CpGtools description for a clock."""
    if name == "EPM":
        return helpdoc.epm_help
    return clockinfo(CLOCK_MODELS[name])


def _add_common_arguments(parser, *, include_log=True):
    """Add arguments shared by most clock subcommands."""
    parser.add_argument(
        "input", type=str, metavar="Input_file", help=helpdoc.input_help
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="out_prefix", default=None,
        help=helpdoc.output_help,
    )
    parser.add_argument(
        "-p", "--percent", type=float, default=0.2, help=helpdoc.na_help
    )
    parser.add_argument(
        "-d", "--delimiter", type=str, default=None, help=helpdoc.del_help
    )
    parser.add_argument(
        "-f", "--format", type=str, choices=("pdf", "png"), default="pdf",
        help=helpdoc.format_help,
    )
    parser.add_argument(
        "-m", "--metadata", type=str, metavar="meta_file", default=None,
        help=helpdoc.meta_help,
    )
    if include_log:
        parser.add_argument(
            "-l", "--log", type=str, metavar="log_file", default=None,
            help=helpdoc.log_help,
        )
    parser.add_argument(
        "--impute", type=int, choices=range(-1, 12), default=11,
        help=helpdoc.imputation_help,
    )
    parser.add_argument(
        "-r", "--ref", type=str, metavar="ref_file", default=None,
        help=helpdoc.ext_ref_help,
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="If set, over-write existing output files.",
    )
    parser.add_argument(
        "--debug", action="store_true", help=helpdoc.debug_help
    )


def _add_mouse_arguments(parser):
    _add_common_arguments(parser)
    parser.add_argument(
        "-g", "--genome", type=str, choices=("mm10", "mm39"),
        default="mm10",
        help=(
            "Reference genome used for mouse RRBS/WGBS read alignment. "
            "Must be 'mm10' or 'mm39'."
        ),
    )


def _add_mammalian_species_arguments(parser):
    _add_common_arguments(parser)
    parser.add_argument(
        "-s", "--species", type=str, choices=("human", "mouse"),
        default="human",
        help="Mammalian species. Currently supports 'human' or 'mouse'.",
    )


def _add_epm_arguments(parser):
    """Add the EPM-specific interface without changing its current options."""
    parser.add_argument(
        "input", type=str, metavar="Input_file", help=helpdoc.input_help
    )
    parser.add_argument(
        "meta", type=str, metavar="meta_file", help=helpdoc.epm_meta_help
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="out_prefix", default=None,
        help=helpdoc.epm_output_help,
    )
    parser.add_argument(
        "-p", "--pcc", type=float, default=0.85,
        help=(
            "Threshold of absolute Pearson correlation coefficient between "
            "chronological age and beta values."
        ),
    )
    parser.add_argument(
        "-n", "--niter", type=int, default=100,
        help="Number of expectation-maximization iterations.",
    )
    parser.add_argument(
        "-k", "--kfold", type=int, default=10,
        help="Number of cross-validation folds.",
    )
    parser.add_argument(
        "-e", "--etol", type=float, default=1e-5,
        help="Error tolerance during model fitting.",
    )
    parser.add_argument(
        "-d", "--delimiter", type=str, default=None, help=helpdoc.del_help
    )
    parser.add_argument(
        "-f", "--format", type=str, choices=("pdf", "png"), default="pdf",
        help=helpdoc.format_help,
    )
    parser.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=helpdoc.log_help,
    )
    parser.add_argument(
        "--impute", type=int, choices=range(-1, 12), default=11,
        help=helpdoc.imputation_help,
    )
    parser.add_argument(
        "-r", "--ref", type=str, metavar="ref_file", default=None,
        help=helpdoc.ext_ref_help,
    )
    parser.add_argument(
        "--debug", action="store_true", help=helpdoc.debug_help
    )


def build_parser():
    """Build and return the epical argument parser."""
    parser = argparse.ArgumentParser(
        description=helpdoc.general_help,
        epilog="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"epical {__version__}"
    )

    subparsers = parser.add_subparsers(
        dest="command", metavar="COMMAND", help="Sub-command description:"
    )

    # Most clocks share the same CLI.
    for name in STANDARD_CLOCKS:
        subparser = subparsers.add_parser(name, help=_clock_help(name))
        _add_common_arguments(subparser)

    # DunedinPACE intentionally has no --log option in the existing CLI.
    dunedin = subparsers.add_parser(
        "DunedinPACE", help=_clock_help("DunedinPACE")
    )
    _add_common_arguments(dunedin, include_log=False)

    epm = subparsers.add_parser("EPM", help=_clock_help("EPM"))
    _add_epm_arguments(epm)

    for name in MOUSE_CLOCKS:
        subparser = subparsers.add_parser(name, help=_clock_help(name))
        _add_mouse_arguments(subparser)

    for name in MAMMALIAN_CLOCKS:
        subparser = subparsers.add_parser(name, help=_clock_help(name))
        _add_mammalian_species_arguments(subparser)

    return parser


def _common_clock_kwargs(args, command):
    """Translate parsed CLI arguments to the methylclocks API."""
    return {
        "beta_file": args.input,
        "outfile": args.output,
        "metafile": args.metadata,
        "delimiter": args.delimiter,
        "cname": command,
        "ff": args.format,
        "na_percent": args.percent,
        "ovr": args.overwrite,
        "imputation_method": args.impute,
        "ext_file": args.ref,
    }


def _configure_logging(args, *, include_logfile=True):
    """Configure CpGtools logging for a parsed command."""
    from dmc.utils import config_log

    if include_logfile:
        config_log(switch=args.debug, logfile=getattr(args, "log", None))
    else:
        config_log(switch=args.debug)


def _run_horvath(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_horvath(
        **_common_clock_kwargs(args, args.command), adult_age=20
    )


def _run_ped_wu(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_horvath(
        **_common_clock_kwargs(args, args.command), adult_age=48
    )


def _run_levine_hannum(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_levine_hannum(
        **_common_clock_kwargs(args, args.command)
    )


def _run_zhang(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_blup_en(**_common_clock_kwargs(args, args.command))


def _run_gestational(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_GA(**_common_clock_kwargs(args, args.command))


def _run_altum_age(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.altum_age(**_common_clock_kwargs(args, args.command))


def _run_epm(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_epm(
        beta_file=args.input,
        metafile=args.meta,
        outfile=args.output,
        delimiter=args.delimiter,
        imputation_method=args.impute,
        ext_file=args.ref,
        pcc_cut=args.pcc,
        iter_n=args.niter,
        error_tol=args.etol,
        cv_folds=args.kfold,
        frmt=args.format,
        cname=args.command,
    )


def _run_general(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_general(**_common_clock_kwargs(args, args.command))


def _run_mammalian(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_mammalian(**_common_clock_kwargs(args, args.command))


def _run_mammalian_species(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_mammalian(
        **_common_clock_kwargs(args, args.command), species=args.species
    )


def _run_mouse(args):
    from dmc import methylclocks

    _configure_logging(args)
    methylclocks.clock_mouse(
        **_common_clock_kwargs(args, args.command), genome=args.genome
    )


def _run_dunedinpace(args):
    from dmc import DunedinPACE

    _configure_logging(args, include_logfile=False)
    kwargs = _common_clock_kwargs(args, args.command)
    kwargs.pop("cname")
    DunedinPACE.DunedinPACE_clock(**kwargs)


def _build_command_handlers():
    """Build the command-to-handler registry.

    Keeping this mapping in one place makes it obvious which implementation
    backs each CLI command and avoids a long conditional dispatch chain.
    """
    handlers = {
        "Ped_Wu": _run_ped_wu,
        "AltumAge": _run_altum_age,
        "EPM": _run_epm,
        "mammClock1": _run_mammalian,
        "DunedinPACE": _run_dunedinpace,
    }

    for command in HORVATH_CLOCKS:
        handlers[command] = _run_horvath
    for command in LEVINE_HANNUM_CLOCKS:
        handlers[command] = _run_levine_hannum
    for command in ZHANG_CLOCKS:
        handlers[command] = _run_zhang
    for command in GESTATIONAL_CLOCKS:
        handlers[command] = _run_gestational
    for command in GENERAL_CLOCKS:
        handlers[command] = _run_general
    for command in MAMMALIAN_CLOCKS:
        handlers[command] = _run_mammalian_species
    for command in MOUSE_CLOCKS:
        handlers[command] = _run_mouse

    return handlers


COMMAND_HANDLERS = _build_command_handlers()


def run_command(args):
    """Dispatch a parsed command through the command-handler registry."""
    try:
        handler = COMMAND_HANDLERS[args.command]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported epical command: {args.command}"
        ) from exc

    handler(args)
    return 0

def main(argv=None):
    """Entry point for the ``epical`` command."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help(sys.stderr)
        return 0

    return run_command(args)


if __name__ == "__main__":
    raise SystemExit(main())

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


# BioLearn models exposed as dedicated epical subcommands.  These are kept
# separate from CpGtools-native clocks so overlapping biological concepts do
# not change the implementation used by existing commands.
BIOLEARN_CLOCKS = (
    "YingCausAge",
    "YingDamAge",
    "YingAdaptAge",
    "DunedinPoAm38",
    "GrimAgeV1",
    "GrimAgeV2",
    "VidalBralo",
    "AlcoholMcCartney",
    "HRSInCHPhenoAge",
    "EpiTOC1",
    "EpiTOC2",
    "SmokingMcCartney",
    "DownSyndrome",
    "StocZ",
    "StocP",
    "StocH",
    "BMI_McCartney",
    "EducationMcCartney",
    "TotalCholesterolMcCartney",
    "HDLCholesterolMcCartney",
    "LDLCholesterolMcCartney",
    "BodyFatMcCartney",
    "BMI_Reed",
    "ProstateCancerKirby",
    "HepatoXu",
    "CVD_Westerman",
    "AD_Bahado-Singh",
    "DepressionBarbu",
    "Bocklandt",
)


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
        "--overwrite", action="store_true",
        help="If set, over-write existing output files.",
    )
    parser.add_argument(
        "--debug", action="store_true", help=helpdoc.debug_help
    )




def _add_biolearn_arguments(parser):
    """Add arguments used by BioLearn-backed subcommands."""
    parser.add_argument(
        "input",
        type=str,
        metavar="Input_file",
        help="Methylation beta-value table (CpGs as rows, samples as columns).",
    )
    parser.add_argument(
        "-m",
        "--metadata",
        type=str,
        metavar="meta_file",
        required=True,
        help=(
            "Sample metadata table. Rows must be samples matching the beta "
            "matrix columns. 'Age'/'Sex' columns are normalized to 'age'/'sex'."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        metavar="out_prefix",
        required=True,
        help=(
            "Output prefix. Results are written to "
            "PREFIX_<BioLearn_model>.tsv."
        ),
    )


def _add_mouse_arguments(parser):
    _add_common_arguments(parser)
    parser.add_argument(
        "-g", "--genome", type=str, choices=("mm10", "mm39"),
        default="mm10",
        help=(
            "Reference genome used for mouse RRBS/WGBS read alignment. "
            "Must be 'mm10' or 'mm39'."),)


def _add_mammalian_species_arguments(parser):
    _add_common_arguments(parser)
    parser.add_argument(
        "-s", "--species", type=str, choices=("human", "mouse"),
        default="human",
        help="Mammalian species. Currently supports 'human' or 'mouse'.",)


def _add_epm_arguments(parser):
    """Add arguments specific to the EPM subcommand."""
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
        "--debug", action="store_true", help=helpdoc.debug_help
    )


def _add_gp_age_arguments(parser):
    """Add arguments specific to the GP-Age subcommand."""
    parser.add_argument(
        "input", type=str, metavar="Input_file", help=helpdoc.input_help
    )
    parser.add_argument(
        "-m", "--model", choices=("10", "30", "71", "a", "b", "c"),
        default="30",
        help="GP-Age model to use.",
    )
    parser.add_argument(
        "-a", "--age", type=str, metavar="age_file", default=None,
        help="Optional file containing chronological ages.",
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="out_prefix", default=None,
        help=(
            "Output prefix. Predictions are written to "
            "PREFIX.GP_age.txt and statistics to PREFIX.GP_age.stats.txt. "
            "If omitted, predictions are written to stdout."
        ),
    )
    parser.add_argument(
        "--model-dir", type=str, metavar="DIR", default=None,
        help=(
            "Directory containing GP-Age model and CpG files. "
            "If omitted, GP-Age uses its bundled model-data directory."
        ),
    )
    parser.add_argument(
        "-l", "--log", type=str, metavar="log_file", default=None,
        help=helpdoc.log_help,
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

    gp_age = subparsers.add_parser(
        "GP-Age",
        help="Gaussian Process-based epigenetic age predictor.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    _add_gp_age_arguments(gp_age)

    epm = subparsers.add_parser("EPM", help=_clock_help("EPM"))
    _add_epm_arguments(epm)

    for name in MOUSE_CLOCKS:
        subparser = subparsers.add_parser(name, help=_clock_help(name))
        _add_mouse_arguments(subparser)

    for name in MAMMALIAN_CLOCKS:
        subparser = subparsers.add_parser(name, help=_clock_help(name))
        _add_mammalian_species_arguments(subparser)


    # Selected BioLearn models.  They intentionally use their BioLearn model
    # names as subcommands, even when CpGtools provides a related clock.
    for name in BIOLEARN_CLOCKS:
        subparser = subparsers.add_parser(
            name,
            help=helpdoc.BIOLEARN_CLOCK_HELP[name],
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        _add_biolearn_arguments(subparser)

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




def _build_biolearn_geodata(beta_path, meta_path):
    """Build a BioLearn GeoData object from beta and metadata tables."""
    import pandas as pd
    from biolearn.data_library import GeoData

    beta_df = pd.read_csv(beta_path, sep=None, index_col=0, engine="python")
    meta_df = pd.read_csv(meta_path, sep=None, index_col=0, engine="python")
    meta_df.rename(columns={"Sex": "sex", "Age": "age"}, inplace=True)

    common_samples = meta_df.index.intersection(beta_df.columns)
    if len(common_samples) == 0:
        raise ValueError(
            "No matching sample IDs found between metadata and beta values."
        )

    meta_df = meta_df.loc[common_samples]
    beta_df = beta_df.loc[:, common_samples]
    return GeoData(metadata=meta_df, dnam=beta_df)


def _run_biolearn(args):
    """Run one of the selected BioLearn models."""
    import pandas as pd
    from biolearn.model_gallery import ModelGallery

    command = args.command
    print(
        f"Reading '{args.input}' and '{args.metadata}' ...",
        file=sys.stderr,
    )
    geodata = _build_biolearn_geodata(args.input, args.metadata)

    print(f"Calculating '{command}' ...", file=sys.stderr)
    output = ModelGallery().get(command).predict(geodata)

    # BioLearn models may return a Series, a one-column DataFrame, or a
    # multi-column DataFrame (e.g. GrimAge and deconvolution models).
    if isinstance(output, pd.Series):
        output = output.to_frame(name=command)
    else:
        output = output.copy()

    output.index.name = "sampleID"

    if command in {"GrimAgeV1", "GrimAgeV2"}:
        renamed = [str(col).replace("DNAm", f"{command}.") for col in output.columns]
        if len(renamed) >= 2:
            renamed[-2] = f"{command}.Prediction"
            renamed[-1] = f"{command}.AgeAcce"
        output.columns = renamed
    elif len(output.columns) == 1:
        output.columns = [command]

    outfile = f"{args.output}_{command}.tsv"
    output.to_csv(outfile, index=True, sep="\t")
    print(f"Wrote '{outfile}'.", file=sys.stderr)



def _run_gp_age(args):
    """Run GP-Age through the CpGtools command-line interface."""
    from dmc import GP_Age

    _configure_logging(args)

    GP_Age.run_gp_age(
        methylation_file=args.input,
        model=args.model,
        age_file=args.age,
        output_prefix=args.output,
        model_dir=args.model_dir,
    )


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
        "GP-Age": _run_gp_age,
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
    for command in BIOLEARN_CLOCKS:
        handlers[command] = _run_biolearn

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

#!/usr/bin/env python3
"""DNA methylation missing-value imputation utilities.

The input matrix is expected to have CpGs in rows and samples in columns.
This refactoring preserves the original imputation algorithms and matrix
orientation while improving validation, structure, logging, and error handling.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import json
from collections.abc import Mapping, Sequence

from cpgmodule._version import __version__
from impute_module.imputelib import fill_statistic,fill_random,fill_ref
from impute_module.imputelib import report_missing_values, toy_df
from impute_module.imputelib import fill_moving_window, fill_knn, fill_buck
from impute_module.imputelib import fill_rf,fill_soft_impute,insert_na
from impute_module.blockimpute import fill_morel
from impute_module.bootstrap import fill_GNN


__author__ = "Liguo Wang"
__license__ = "MIT"

LOGGER = logging.getLogger("nafiller")

def load_morel_groups(
    path: str | Path | None,) -> Mapping[str, Sequence[str]] | None:
    """
    Load the two MOREL sample groups from a JSON file.

    Expected format
    ---------------
    {
        "group_1": ["sample1", "sample2"],
        "group_2": ["sample3", "sample4"]
    }
    """
    if path is None:
        return None

    group_path = Path(path)

    try:
        with group_path.open(encoding="utf-8") as handle:
            groups = json.load(handle)
    except OSError as exc:
        raise ValueError(
            f"Unable to read group file {group_path!s}: {exc}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON in group file {group_path!s}: {exc}"
        ) from exc

    if not isinstance(groups, dict):
        raise ValueError(
            "The group JSON must contain an object mapping group names "
            "to sample-name lists."
        )

    return groups

def configure_logging(debug: bool = False) -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if debug else logging.INFO,)

def read_df(infile: str | Path) -> pd.DataFrame:
    """
    Read a DNA methylation matrix with CpGs in rows and samples in columns.
    """
    path = Path(infile)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    LOGGER.info('Reading input file "%s" as a data frame.', path)
    df = pd.read_csv(path, index_col=0, sep=None, engine="python", compression="infer",)
    if df.empty:
        raise ValueError(f"Input matrix is empty: {path}")
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def write_df(df: pd.DataFrame, outfile: str | Path, decimal: int = 5,
    overwrite: bool = False, na_repr: str = "NA",) -> None:
    """Round and write a matrix as a tab-delimited file."""
    if decimal < 0:
        raise ValueError("decimal must be a non-negative integer")

    path = Path(outfile)
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output file already exists: {path}. Use --overwrite to replace it."
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    rounded = df.round(decimal)
    rounded.to_csv(path, sep="\t", na_rep=na_repr)
    LOGGER.info(
        'File "%s" written successfully (%d missing values).',path,
        int(rounded.isna().sum().sum()),)

def add_io_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("input", metavar="input_file",
        help="Input matrix; CpGs in rows, samples in columns.")
    parser.add_argument("output", metavar="output_file",
        help="Output tab-delimited matrix.")
    parser.add_argument("--decimal",type=int,default=5,
        help="Number of decimal places in the output.",)
    parser.add_argument("--overwrite", action="store_true",
        help="Overwrite an existing output file.",)


def add_truth_argument(parser: argparse.ArgumentParser) -> None:
    """Add an optional truth-matrix argument to an imputation command."""
    parser.add_argument(
        "--truth",
        metavar="FILE",
        help=(
            "Optional complete or partially complete truth matrix used to "
            "evaluate imputed values. Metrics are calculated only at cells "
            "that were missing in the original input and are observed in "
            "both the truth and imputed matrices."
        ),
    )


def evaluate_imputation(
    input_df: pd.DataFrame,
    imputed_df: pd.DataFrame,
    truth_file: str | Path,
) -> None:
    """Report MAE, RMSE, and coefficient of determination (R2).

    Evaluation is restricted to cells that were missing in ``input_df``.
    Cells that are missing in either the truth matrix or the imputed matrix
    are excluded, which allows partially imputed output matrices.
    """
    truth_path = Path(truth_file)
    if not truth_path.is_file():
        raise FileNotFoundError(f"Truth file not found: {truth_path}")

    LOGGER.info('Reading truth file "%s" as a data frame.', truth_path)
    truth_df = pd.read_csv(
        truth_path,
        index_col=0,
        sep=None,
        engine="python",
        compression="infer",
    )

    if truth_df.empty:
        raise ValueError(f"Truth matrix is empty: {truth_path}")

    truth_df.index = truth_df.index.astype(str)
    truth_df.columns = truth_df.columns.astype(str)
    truth_df = truth_df.apply(pd.to_numeric, errors="coerce")

    missing_rows = [name for name in input_df.index if name not in truth_df.index]
    missing_cols = [name for name in input_df.columns if name not in truth_df.columns]
    if missing_rows or missing_cols:
        details = []
        if missing_rows:
            preview = missing_rows[:10]
            details.append(
                f"{len(missing_rows)} CpG IDs missing from truth, e.g. {preview}"
            )
        if missing_cols:
            preview = missing_cols[:10]
            details.append(
                f"{len(missing_cols)} sample IDs missing from truth, e.g. {preview}"
            )
        raise ValueError(
            "Truth matrix does not contain all input rows and columns: "
            + "; ".join(details)
        )

    missing_output_rows = [
        name for name in input_df.index if name not in imputed_df.index
    ]
    missing_output_cols = [
        name for name in input_df.columns if name not in imputed_df.columns
    ]
    if missing_output_rows or missing_output_cols:
        raise ValueError(
            "Imputed matrix does not retain all rows and columns from the "
            "original input."
        )

    truth_aligned = truth_df.reindex(
        index=input_df.index,
        columns=input_df.columns,
    )
    imputed_aligned = imputed_df.reindex(
        index=input_df.index,
        columns=input_df.columns,
    )

    original_missing = input_df.isna().to_numpy()
    truth_values = truth_aligned.to_numpy(dtype=float)
    imputed_values = imputed_aligned.to_numpy(dtype=float)

    truth_available = np.isfinite(truth_values)
    imputed_available = np.isfinite(imputed_values)
    evaluable = original_missing & truth_available & imputed_available

    total_original_missing = int(original_missing.sum())
    truth_missing_count = int((original_missing & ~truth_available).sum())
    still_missing_count = int(
        (original_missing & truth_available & ~imputed_available).sum()
    )
    evaluated_count = int(evaluable.sum())

    LOGGER.info(
        "Truth evaluation: %d of %d originally missing values are evaluable.",
        evaluated_count,
        total_original_missing,
    )

    if truth_missing_count:
        LOGGER.info(
            "Excluded %d positions because the truth matrix is missing.",
            truth_missing_count,
        )
    if still_missing_count:
        LOGGER.info(
            "Excluded %d positions because the imputed output is still missing.",
            still_missing_count,
        )

    if evaluated_count == 0:
        LOGGER.warning(
            "No originally missing positions have values in both the truth "
            "and imputed matrices; metrics were not calculated."
        )
        return

    y_true = truth_values[evaluable]
    y_pred = imputed_values[evaluable]
    residual = y_pred - y_true

    mae = float(np.mean(np.abs(residual)))
    rmse = float(np.sqrt(np.mean(np.square(residual))))

    ss_res = float(np.sum(np.square(residual)))
    centered = y_true - float(np.mean(y_true))
    ss_tot = float(np.sum(np.square(centered)))

    if evaluated_count < 2 or ss_tot == 0.0:
        r2 = float("nan")
        r2_text = "NA (fewer than two distinct truth values)"
    else:
        r2 = 1.0 - (ss_res / ss_tot)
        r2_text = f"{r2:.6f}"

    LOGGER.info(
        "Imputation accuracy on %d values: MAE=%.6f, RMSE=%.6f, R2=%s.",
        evaluated_count,
        mae,
        rmse,
        r2_text,
    )

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=(
            "Impute missing values in DNA methylation matrices. CpGs must be "
            "in rows and samples must be in columns."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    parser.add_argument("-v", "--version", action="version",
        version=f"nafiller {__version__}",)
   
    parser.add_argument("--debug", action="store_true", 
        help="Enable debug logging.")

    subparsers = parser.add_subparsers(dest="command", required=True)
    
    ###################
    ## Toy df
    ###################
    p = subparsers.add_parser("toy",
        help="Generate a toy data frame with missing values.",
        description=(
            "Generate a random matrix with user-specified "
            "dimensions and randomly inserted missing values (NA/NaN)."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    p.add_argument("output", metavar="output_file",
        help="Output file for the generated data frame.",)
    p.add_argument("--rows", type=int, default=20, metavar="INT",
        help="Number of rows (CpGs).",)
    p.add_argument("--cols", type=int, default=5, metavar="INT",
        help="Number of columns (samples).",)
    p.add_argument( "--missingness", type=float, default=0.2, metavar="FLOAT",
        help=(
            "Amount of missing data to introduce. Values between 0 and 1 "
            "specify the fraction of matrix entries to replace with NA; "
            "values greater than or equal to 1 specify the approximate "
            "number of missing values."),)
    p.add_argument("--min", dest="min_val", type=float, default=0.0,
        metavar="FLOAT",
        help="Minimum random value.",)
    p.add_argument("--max", dest="max_val", type=float, default=1.0,
        metavar="FLOAT",
        help="Maximum random value.",)
    p.add_argument("--seed", type=int, default=1234, metavar="INT",
        help="Random seed.",)
    p.add_argument("--prefix", default="sample", metavar="STR",
        help="Prefix used for sample (column) names.",)
    p.add_argument("--decimal", type=int, default=5, metavar="INT",
        help="Number of decimal places in the output matrix.",)
    p.add_argument("--na-repr", default="NA", metavar="STR",
        help="String used to represent missing values in the output file.",)


    ############################
    ## Insert missing values
    ############################
    p = subparsers.add_parser(
        "insertna",
        help="Randomly insert missing values into a matrix.",
        description=(
            "Randomly replace observed matrix entries with missing values until "
            "the requested total number of missing values is reached. Existing "
            "missing values are preserved."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    
    p.add_argument(
        "target_missing",
        type=int,
        metavar="INT",
        help=(
            "Target total number of missing values in the output matrix. "
            "For example, if the input already contains 20 missing values and "
            "TARGET_MISSING is 100, 80 additional missing values are inserted."
        ),)
    
    p.add_argument(
        "--seed",
        dest="random_state",
        type=int,
        default=123,
        metavar="INT",
        help="Random seed used to select positions for new missing values.",)

    p.add_argument(
        "--na-repr",
        default="NA",
        metavar="STR",
        help="String used to represent missing values in the output file.",
    )

    ###################
    ## Remove NAs
    ###################
    p = subparsers.add_parser(
        "dropna",
        help="Remove rows or columns containing missing values.",
        description=(
            "Remove rows or columns that contain one or more missing values "
            "(NA/NaN) from input file."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    
    p.add_argument(
        "--axis",
        choices=("rows", "columns"),
        default="rows",
        help=("Remove rows (CpGs) or columns (samples) containing "
            "missing values."),)

    ###################
    ## count NAs
    ###################
    p = subparsers.add_parser(
        "countna",
        help="Count missing values per row (or per column).",
        description=(
            "Count missing values (NA/NaN) in a DNA methylation matrix and "
            "generate two reports: one containing the number of missing values "
            "for each CpG (row) and another containing the number of missing "
            "values for each sample (column)."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    p.add_argument("input", metavar="input_file",
        help="Input DNA methylation matrix.",)
    
    p.add_argument("--row-report", dest="row_na", default="missing_by_cpg.tsv",
        metavar="FILE",
        help=("Output file containing the number of missing values for each "
            "CpG (row)."),)
    
    p.add_argument("--column-report", dest="column_na", 
        default="missing_by_sample.tsv",
        metavar="FILE",
        help=("Output file containing the number of missing values for each "
            "sample (column)."),)

    ###########################
    ## Fill NAs with a constant
    ###########################
    p = subparsers.add_parser(
        "constant",
        help="Replace missing values with a fixed value.",
        description=(
            "Replace all missing values (NA/NaN) with a user-specified "
            "constant."
        ),formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    add_io_arguments(p)
    add_truth_argument(p)
    p.add_argument(
        "x",
        nargs="?",
        type=float,
        default=0.0,
        metavar="VALUE",
        help=("Value used to replace all missing values. "
            "Default: %(default)s."),)

    ###########################
    ## Fill NAs with statistics
    ###########################
    for name, help_text in (
        ("mean", "Impute missing values using row or column means.",),
        ("median", "Impute missing values using row or column medians.",),
        ("min", "Impute missing values using row or column minima.",),
        ("max", "Impute missing values using row or column maxima.",),
        ("rand", 
         "Impute using randomly selected values from the same row or column.",),
        ):
        p = subparsers.add_parser(
            name,
            help=help_text,
            description=help_text,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
        add_io_arguments(p)
        add_truth_argument(p)
    
        p.add_argument(
            "--axis",
            choices=("index", "columns"),
            default="index",
            metavar="{index,columns}",
            help=(
                "Dimension used to calculate replacement values. "
                "'index' = operate on rows (CpGs); "
                "'columns' = operate on columns (samples)."),)

    ##########################################
    ## Fill NAs with external reference (KNN)
    ##########################################
    p = subparsers.add_parser(
        "refknn",
        help="Impute missing values using nearest neighbors from an external reference matrix.",
        description=(
            "Impute missing values with a K-nearest-neighbors algorithm that "
            "searches for similar rows or columns in an external, complete "
            "reference matrix. Missing values are replaced with the mean of "
            "the corresponding values from the selected neighbors."
        ),formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    add_io_arguments(p)
    add_truth_argument(p)
    
    p.add_argument( "-r", "--reference", dest="reference", required=True,
        metavar="FILE", help=(
            "External reference matrix used for nearest-neighbor searches. "
            "The reference matrix must not contain missing values."),)
    
    p.add_argument("--search_axis", choices=("index", "columns"),
        default="columns", metavar="{index,columns}", help=(
            "Dimension in which nearest neighbors are searched. "
            "'columns' searches for similar columns (samples); the reference "
            "must contain the input CpGs. "
            "'index' searches for similar rows (CpGs); the reference must "
            "contain the input samples."),)
    
    p.add_argument("-k", "--neighbors", type=int, default=3,  metavar="INT",
        help="Number of nearest neighbors used for imputation.",)

    ###################
    ## Movign window
    ###################
    p = subparsers.add_parser("mw",
        help="Impute missing values using a moving window.",
        description=(
            "Impute missing values using neighboring values within a moving "
            "window. The window is applied independently along rows (CpGs) or "
            "columns (samples), and neighboring values are summarized using "
            "the mean or median."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    add_io_arguments(p)
    add_truth_argument(p)
    p.add_argument("--axis", choices=("index", "columns"), default="index",
        metavar="{index,columns}", help=(
            "Dimension along which the moving window is applied. "
            "'index' = apply a horizontal moving window to each row; "
            "'columns' = apply a vertical moving window to each column."),)
    p.add_argument("--naindex", type=int, choices=(0, -1), default=None,
        metavar="{0,-1}", help=(
            "Position of the missing value within the moving window. "
            "0 = use values to the right only; "
            "-1 = use values to the left only. "
            "If omitted, center the missing value within the window whenever "
            "possible."),)
    p.add_argument("--wsize", type=int, default=5, metavar="INT",
        help=(
            "Size of the moving window, including the position of the missing "
            "value."),)
    p.add_argument("--errors",
        choices=("raise", "coerce", "ignore"),
        default="coerce",
        metavar="{raise,coerce,ignore}",
        help=("How to handle windows that extend beyond the matrix boundary. "
            "'raise' = raise an exception; "
            "'coerce' = retry with a centered window; "
            "'ignore' = leave the missing value unchanged."),)
    p.add_argument("--func", choices=("mean", "median"), default="mean",
        metavar="{mean,median}",
        help="Summary statistic used to estimate missing values.",)
    
    ###################
    ## KNN
    ###################    
    p = subparsers.add_parser(
        "knn",
        help="Impute missing values using scikit-learn's KNNImputer.",
        description=(
            "Impute missing values using K-nearest neighbors within the input "
            "matrix. Neighbors can be searched among columns (samples) or rows "
            "(CpGs). Distances are calculated using scikit-learn's "
            "NaN-aware Euclidean distance, so initial mean imputation is not "
            "required."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    add_truth_argument(p)
    p.add_argument(
        "--axis",
        choices=("index", "columns"),
        default="columns",
        metavar="{index,columns}",
        help=(
            "Dimension in which nearest neighbors are searched. "
            "'columns' = search neighboring samples; "
            "'index' = search neighboring CpGs."),)
    p.add_argument(
        "-k",
        "--neighbors",
        type=int,
        default=3,
        metavar="INT",
        help="Number of nearest neighbors used for imputation.",)
    p.add_argument(
        "--weights",
        choices=("uniform", "distance"),
        default="distance",
        metavar="{uniform,distance}",
        help=(
            "Weighting method used when averaging neighbor values. "
            "'uniform' gives every neighbor equal weight; "
            "'distance' gives closer neighbors greater weight."),)
    p.add_argument(
        "--keep-empty-features",
        action="store_true",
        help=(
            "Retain features that contain only missing values. Their imputed "
            "values are set to zero by KNNImputer."),)

    ###################
    ## Buck interative
    ################### 
    p = subparsers.add_parser(
        "buck",
        help="Impute missing values using iterative Buck regression.",
        description=(
            "Impute missing values using the iterative regression method of "
            "Buck (1960). Missing values are initially replaced by row "
            "means and then iteratively updated using linear regression until "
            "the estimates converge."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    add_truth_argument(p)
    p.add_argument("--eps", type=float, default=1e-3, metavar="FLOAT",
        help=(
            "Convergence threshold. Iteration stops when the maximum relative "
            "change in the imputed values is less than this value."),)
    p.add_argument("--max-iter", type=int, default=100, metavar="INT",
        help="Maximum number of iterations.",)


    ############################
    ## Random-forest iterative
    ############################
    p = subparsers.add_parser(
        "rf",
        help="Impute missing values using iterative random-forest regression.",
        description=(
            "Impute missing values using an iterative random-forest procedure. "
            "Missing values are initialized using row means. Each incomplete "
            "column is then predicted from the remaining columns using "
            "RandomForestRegressor, and the process is repeated until convergence "
            "or until the maximum number of iterations is reached. "
            "For beta-value matrices, specify --min-value 0.0 and --max-value "
            "1.0. For M-value matrices, leave both options unspecified."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    add_io_arguments(p)
    add_truth_argument(p)
    p.add_argument("--trees", dest="n_estimators", type=int, default=100,
        metavar="INT",
        help="Number of trees fitted for each incomplete column.",)
    p.add_argument("--max-iter", type=int, default=10, metavar="INT",
        help="Maximum number of complete random-forest imputation rounds.",)
    p.add_argument("--eps", type=float, default=1e-3, metavar="FLOAT",
        help=(
            "Convergence threshold. Iteration stops when the maximum relative "
            "change among imputed values is smaller than this value."),)
    p.add_argument("--max-depth", type=int, default=None, metavar="INT",
        help=(
            "Maximum tree depth. If omitted, trees expand until another stopping "
            "criterion is reached."),)
    p.add_argument("--min-samples-leaf", type=int, default=1, metavar="INT",
        help="Minimum number of observations required in a terminal leaf.",)
    p.add_argument("--seed", dest="random_state", type=int, default=1234,
        metavar="INT",
        help="Random seed used for reproducible forest fitting.",)
    p.add_argument("--jobs", dest="n_jobs", type=int, default=-1,metavar="INT",
        help=(
            "Number of parallel jobs used to fit each forest. "
            "-1 uses all available processors."),)
    p.add_argument("--min-value", type=float, default=0.0, metavar="FLOAT",
        help="Minimum permitted imputed value.",)
    p.add_argument("--max-value", type=float, default=1.0, metavar="FLOAT",
        help="Maximum permitted imputed value.",)


    ###################
    ## SoftImpute
    ###################
    p = subparsers.add_parser(
        "softimpute",
        help="Impute missing values using low-rank SoftImpute matrix completion.",
        description=(
            "Impute missing values by repeatedly applying soft-thresholded "
            "singular-value decomposition. SoftImpute assumes that the DNA "
            "methylation matrix can be approximated by a low-rank matrix. "
            "Observed values are preserved while missing values are updated "
            "until convergence or the maximum iteration count is reached. "
            "For beta-value matrices, specify --min-value 0.0 and --max-value "
            "1.0. For M-value matrices, leave both options unspecified."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    add_truth_argument(p)
    
    p.add_argument("--shrinkage", dest="shrinkage_value", type=float,
        default=None,
        metavar="FLOAT",
        help=(
            "Nonnegative value subtracted from each singular value. "
            "If omitted, the value is estimated from the largest singular "
            "value of the initialized matrix."),)
    p.add_argument("--threshold", dest="convergence_threshold", type=float,
        default=1e-3,
        metavar="FLOAT",
        help=(
            "Convergence threshold for the relative change among the imputed "
            "values."),)
    p.add_argument("--max-iter",dest="max_iters",type=int,default=100,
        metavar="INT",
        help="Maximum number of SoftImpute iterations.",)
    p.add_argument("--max-rank",type=int,default=10, metavar="INT",
        help=(
            "Maximum rank used for truncated randomized SVD. A value between "
            "5 and 20 is a reasonable starting range for large methylation "
            "matrices."),)
    p.add_argument("--power-iter",dest="n_power_iterations",type=int,default=2,
        metavar="INT",
        help=(
            "Number of power iterations used by randomized SVD. Larger values "
            "can improve decomposition accuracy but increase runtime."),)
    p.add_argument("--init-fill",dest="init_fill_method",
        choices=("zero", "mean", "median", "min", "random"),
        default="zero",
        metavar="{zero,mean,median,min,random}",
        help="Method used to initialize missing entries before iteration.",)
    p.add_argument("--min-value",type=float,default=None,metavar="FLOAT",
        help=(
            "Minimum permitted imputed value. The default is appropriate for "
            "methylation beta values."),)
    p.add_argument("--max-value",type=float,default=None,
        metavar="FLOAT",
        help=(
            "Maximum permitted imputed value. The default is appropriate for "
            "methylation beta values."),)
    p.add_argument("--seed",dest="random_state",type=int,default=1234,
        metavar="INT",
        help="Random seed used by randomized SVD.",)
    
    p.add_argument("--verbose-softimpute",action="store_true",
        help="Log SoftImpute progress for every iteration.",)


    ###################
    ## MOREL imputation
    ###################
    p = subparsers.add_parser(
        "morel",
        help="Impute block-wise missing values using the MOREL method.",
        description=(
            "Impute systematic block-wise missing values in a DNA "
            "methylation matrix. Samples are divided into two groups, and "
            "CpGs observed in one group but entirely missing in the other "
            "are predicted using a Random Forest or dense neural network. "
            "Sporadic missing values are first imputed using K-nearest "
            "neighbors. For beta-value matrices, specify --min-value 0.0 and "
            "--max-value 1.0. For M-value matrices, leave both options unspecified."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    add_truth_argument(p)
    
    p.add_argument("--group", metavar="FILE",
        help=(
            "Optional file defining exactly two sample groups. Each line "
            "should contain a group name followed by a comma-separated "
            "list of sample names. If omitted, the two groups are inferred "
            "from sample missingness patterns using K-means clustering."),)
    p.add_argument("--model", choices=("RF", "DNN"), default="RF",
        help=(
            "Secondary model used to predict block-wise missing values: "
            "'RF' for Random Forest or 'DNN' for a dense neural network."),)
    p.add_argument("--knn-neighbors", type=int, default=5, metavar="INT",
        help=(
            "Number of nearest neighbors used to impute sporadic "
            "missing values."),)
    p.add_argument("--knn-weights", choices=("uniform", "distance"),
        default="uniform",
        help="Weighting strategy used by K-nearest-neighbor imputation.",)
    p.add_argument("--n-iter", type=int, default=10, metavar="INT",
        help=(
            "Number of repeated model fits when using Random Forest. "
            "Predictions from all iterations are averaged."),)
    p.add_argument("--seed", dest="random_state", type=int,
        default=100, metavar="INT",
        help="Random seed used for data splitting and model fitting.",)
    p.add_argument("--n-jobs", type=int, default=-1, metavar="INT",
        help=(
            "Number of parallel jobs used by Random Forest. "
            "-1 uses all available processors."),)
    p.add_argument("--train-size", type=float, default=0.75,
        metavar="FLOAT",
        help=(
            "Fraction of complete CpGs used for model training during "
            "validation."),)
    p.add_argument("--n-estimators", type=int, default=100, metavar="INT",
        help="Number of trees in each Random Forest.",)
    p.add_argument("--max-depth", type=int, default=30, metavar="INT",
        help="Maximum depth of each Random Forest tree.",)
    p.add_argument("--min-samples-leaf", type=int, default=1, metavar="INT",
        help=(
            "Minimum number of observations required in a terminal "
            "Random Forest leaf."),)
    p.add_argument("--max-features", type=float, default=1.0, metavar="FLOAT",
        help=(
            "Fraction or number of predictors considered at each "
            "Random Forest split."),)
    p.add_argument("--epochs", type=int, default=100, metavar="INT",
        help=(
            "Maximum number of training epochs when using the dense "
            "neural network."),)
    p.add_argument("--min-value", dest="min_value", type=float, default=0.0,
        metavar="FLOAT",
        help="Lower bound applied to imputed values.",
    )
    
    p.add_argument("--max-value", dest="max_value", type=float,
        default=1.0,
        metavar="FLOAT",
        help="Upper bound applied to imputed values.",
    )
    p.add_argument("--na-repr", default="NA", metavar="STR",
        help=(
            "String used to represent values that remain missing in "
            "the output file."),)

    ###################
    ## Genomic KNN
    ###################
    p = subparsers.add_parser(
        "gnn",
        help="Impute missing values using CpGs that are genomic nearest neighbors.",
        description=(
            "Impute missing DNA methylation values using nearby CpGs in the "
            "genome. Candidate neighbors are selected according to genomic "
            "distance and may optionally be restricted to the same candidate "
            "cis-regulatory element (CRE), such as a promoter, enhancer, or CpG "
            "island."
        ),formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    
    add_io_arguments(p)
    add_truth_argument(p)
    
    p.add_argument(
        "--gfile",
        required=True,
        metavar="FILE",
        help=(
            "CpG genomic annotation file. Must contain at least five columns: "
            "chrom, start, end, cpg_id, and CRE."    ),)
    
    p.add_argument(
        "--cpgfile",
        metavar="FILE",
        help=(
            "Optional file containing additional CpG IDs (one per line). "
            "These CpGs will be appended to the input matrix and imputed when "
            "possible. For example, you might want to add 450K-specific "
            "CpGs to 850K matrix."),)
    
    p.add_argument(
        "--up-dist",
        type=int,
        default=100,
        metavar="INT",
        help=(
            "Maximum genomic distance (bp) for searching upstream CpGs."),)
    
    p.add_argument(
        "--down-dist",
        type=int,
        default=100,
        metavar="INT",
        help=(
            "Maximum genomic distance (bp) for searching downstream CpGs."    ),)
    
    p.add_argument(
        "--up-ncpg",
        type=int,
        default=2,
        metavar="INT",
        help="Maximum number of upstream CpGs used for imputation.",)
    
    p.add_argument(
        "--down-ncpg",
        type=int,
        default=2,
        metavar="INT",
        help="Maximum number of downstream CpGs used for imputation.",)
    
    p.add_argument(
        "--same-cre",
        action="store_true",
        help=(
            "Restrict candidate neighbors to CpGs sharing at least one "
            "candidate regulatory element (CRE) with the target CpG."    ),)
    
    p.add_argument(
        "--method",
        choices=("AA", "WA", "TA"),
        default="TA",
        help=(
            "Method used to combine neighboring beta values: "
            "AA (arithmetic average), "
            "WA (inverse-distance weighted average), or "
            "TA (trimmed average)."    ),)
    
    p.add_argument(
        "--na-repr",
        default="NA",
        metavar="STR",
        help="String used to represent missing values in the output file.",)



    return parser

def run_command(args: argparse.Namespace) -> None:
    command = args.command.lower()
    
    # create a toy df for testing. No input required.
    if command == "toy":
        output = toy_df(
            n_rows=args.rows,
            n_cols=args.cols,
            missingness=args.missingness,
            sample_prefix=args.prefix,
            min_val=args.min_val,
            max_val=args.max_val,
            rand_seed=args.seed,
        )
        output.to_csv(
            args.output,
            sep="\t",
            index=True,
            float_format=f"%.{args.decimal}f",
            na_rep=args.na_repr,
        )
        return
    
    # Read input file.
    input_df = read_df(args.input)
    
    # insert NA, do not need to check the existance of missing. 
    if command == "insertna":
        output_df = insert_na(
            df=input_df,
            target_missing=args.target_missing,
            random_state=args.random_state,)
        write_df(output_df, args.output, args.decimal, args.overwrite,
                 na_repr=args.na_repr)
        return
    
    # Missing-value checks for imputation commands
    na_mask = input_df.isna()
    total_missing = int(na_mask.to_numpy().sum())
    row_missing = int(na_mask.any(axis=1).sum())
    col_missing = int(na_mask.any(axis=0).sum())
    LOGGER.info(
        'Input file "%s" contains %d missing values.',args.input,total_missing)
    LOGGER.info(
        "%d rows contain at least one missing value.",row_missing,)
    LOGGER.info(
        "%d columns contain at least one missing value.",col_missing,)
    if total_missing == 0:
        LOGGER.info("No missing values were found. Nothing to do.")
        sys.exit(0)


    if command == "dropna":        
        LOGGER.info("Removing *%s* containing any missing values.", args.axis)
        axis_name = 0 if args.axis == "rows" else 1
        output_df = input_df.dropna(axis=axis_name, how="any")

    elif command == "countna":        
        LOGGER.info("Save row-wise missing values to \"%s\"", args.row_na)
        LOGGER.info("Save column-wise missing values to \"%s\"", args.column_na)
        report_missing_values(input_df, args.row_na, args.column_na)

    elif command == "constant":
        LOGGER.info("Replacing missing values with %s.", args.x)
        output_df = input_df.fillna(args.x)

    elif command == "mean":
        LOGGER.info("Replacing missing values with %s means.", "column" if args.axis == 'columns' else "row")
        output_df = fill_statistic(input_df, args.axis, "mean")

    elif command == "median":
        LOGGER.info("Replacing missing values with %s medians.", "column" if args.axis == 'columns' else "row")
        output_df = fill_statistic(input_df, args.axis, "median")

    elif command == "min":
        LOGGER.info("Replacing missing values with %s minima.", "column" if args.axis == 'columns' else "row")
        output_df = fill_statistic(input_df, args.axis, "min")

    elif command == "max":
        LOGGER.info("Replacing missing values with %s maxima.", "column" if args.axis == 'columns' else "row")
        output_df = fill_statistic(input_df, args.axis, "max")

    elif command == "rand":
        LOGGER.info("Replacing missing values with random values from the same %s.", "column" if args.axis == 'columns' else "row")
        output_df = fill_random(input_df, args.axis)

    elif command == "refknn":
        LOGGER.info("Replacing missing values from the external reference.")
        # read external ref file as a DataFrame
        refdf = pd.read_csv(args.reference,
            index_col=0,sep=None, engine="python",compression="infer",)
        output_df = fill_ref(df = input_df,
                             ref = refdf, 
                             axis = args.search_axis,
                             k = args.neighbors)
    elif command == "mw":
        LOGGER.info("Impute missing values by applying moving window to each %s.", "column" if args.axis == 'columns' else "row")
        output_df = fill_moving_window(
            df=input_df,
            axis=args.axis,
            naindex=args.naindex,
            errors=args.errors,
            wsize=args.wsize,
            func=args.func,)
    elif command == "knn":
        LOGGER.info("Impute missing values using KNN. Search neighbours from %s.", "column" if args.axis == 'columns' else "row")
        output_df = fill_knn(
            df=input_df,
            axis=args.axis,
            neighbors=args.neighbors,
            weights=args.weights,
            keep_empty_features=args.keep_empty_features,)
    elif command == 'buck':
        LOGGER.info("Impute missing values using Buck's method (Buck, S. F., 1960).")
        output_df = fill_buck(
            df=input_df,
            eps=args.eps,
            max_iter=args.max_iter)
    elif command == 'rf':
        LOGGER.info("Impute missing values using iterative Random Forest regression.")
        output_df = fill_rf(
            df=input_df,
            n_estimators=args.n_estimators,
            max_iter=args.max_iter,
            eps=args.eps,
            max_depth=args.max_depth,
            min_samples_leaf=args.min_samples_leaf,
            random_state=args.random_state,
            n_jobs=args.n_jobs,
            min_value=args.min_value,
            max_value=args.max_value,
            )     
    elif command == 'softimpute':
        LOGGER.info("Impute missing values using SoftImpute.")
        output_df = fill_soft_impute(
            df=input_df,
            shrinkage_value=args.shrinkage_value,
            convergence_threshold=args.convergence_threshold,
            max_rank=args.max_rank,
            max_iters=args.max_iters,
            n_power_iterations=args.n_power_iterations,
            init_fill_method=args.init_fill_method,
            min_value=args.min_value,
            max_value=args.max_value,
            random_state=args.random_state,
            )
    elif command == "morel":
        LOGGER.info(
            "Impute block-wise missing values using MOREL method (%s).",
            args.model,
        )

        group = None
        if args.group is not None:
            group = load_morel_groups(args.group)   # returns dict[str, list[str]]
    
        output_df = fill_morel(
            df=input_df,
            group=group,
            model=args.model,
            decimal=args.decimal,
            knn_neighbors=args.knn_neighbors,
            knn_weights=args.knn_weights,
            n_iter=args.n_iter,
            random_state=args.random_state,
            n_jobs=args.n_jobs,
            train_size=args.train_size,
            n_estimators=args.n_estimators,
            max_depth=args.max_depth,
            min_samples_leaf=args.min_samples_leaf,
            max_features=args.max_features,
            epochs=args.epochs,
            min_value=args.min_value,
            max_value=args.max_value,
        )
    elif command == 'gnn':
        LOGGER.info(
            "Impute missing values using neighboring genomic CpGs.")
        output_df = fill_GNN(
            df = input_df,
            gfile = args.gfile,
            cpgfile = args.cpgfile,
            up_dist = args.up_dist,
            down_dist =  args.down_dist,
            up_ncpg= args.up_ncpg,
            down_ncpg= args.down_ncpg,
            same_CRE = args.same_cre,
            method = args.method,
        )       
        
    else:
        raise ValueError(f"Unknown command: {args.command}")
    
    if command not in ['countna']:
        write_df(output_df, args.output, args.decimal, args.overwrite)

        truth_file = getattr(args, "truth", None)
        if truth_file is not None:
            # write_df rounds values before writing, so evaluate the same values
            # that appear in the output file.
            evaluate_imputation(
                input_df=input_df,
                imputed_df=output_df.round(args.decimal),
                truth_file=truth_file,
            )


def nafiller(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.debug)

    try:
        run_command(args)
    except Exception as exc:
        LOGGER.error("%s", exc)
        if args.debug:
            LOGGER.exception("Detailed traceback")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(nafiller())

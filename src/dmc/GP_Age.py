#!/usr/bin/env python3
"""
Predict biological age using GP-Age.

GP-Age is a Gaussian Process-based epigenetic age predictor using
predefined DNA methylation CpG sites.
"""

import os
import sys
import logging
import argparse

import GPy
import numpy as np
import pandas as pd

from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    mean_squared_error,
    median_absolute_error,
    mean_absolute_error,
)


logger = logging.getLogger(__name__)

MODEL_CHOICES = ("10", "30", "71", "a", "b", "c")


def model_name(model):
    """Convert command-line model name to GP-Age model suffix."""
    if model in ("a", "b", "c"):
        return model
    return "%s_cpgs" % model


def get_model_dir():
    """
    Return directory containing GP-Age model files.

    Assumes model data are installed in a ``data`` directory
    adjacent to this script.
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def load_sites(model, model_dir):
    """Load CpG sites required by a GP-Age model."""
    name = model_name(model)
    site_file = os.path.join(
        model_dir,
        "GP-age_sites_%s.csv" % name,
    )

    if not os.path.isfile(site_file):
        raise FileNotFoundError(
            "Cannot find GP-Age CpG file: %s" % site_file
        )

    df = pd.read_csv(site_file)

    if df.empty:
        raise ValueError(
            "GP-Age CpG file is empty: %s" % site_file
        )

    sites = (
        df.iloc[:, -1]
        .dropna()
        .astype(str)
        .str.strip()
    )

    if len(sites) == 0:
        raise ValueError(
            "No CpG sites were found in: %s" % site_file
        )

    return pd.Index(sites)


def load_model(model, model_dir):
    """Load a pretrained GP-Age model."""
    name = model_name(model)
    model_file = os.path.join(model_dir, "GP-age_model_%s.json.zip" % name)

    if not os.path.isfile(model_file):
        raise FileNotFoundError("Cannot find GP-Age model: %s" % model_file)

    logger.info("Loading GP-Age model: %s", name)
    predictor = GPy.models.GPRegression.load_model(model_file)
    logger.info("GP-Age model loaded")

    return predictor


def read_methylation(infile):
    """
    Read DNA methylation matrix.

    Input format:
        rows    = CpG sites
        columns = samples

    The first column must contain CpG identifiers.
    """
    logger.info("Reading methylation data: %s", infile)

    meth = pd.read_csv(
        infile,
        sep=None,
        engine="python",
        index_col=0,
    )

    if meth.empty:
        raise ValueError("Methylation input file is empty")

    # GP-Age expects samples in rows and CpGs in columns.
    meth = meth.T

    logger.info(
        "Loaded %d samples and %d CpG sites",
        meth.shape[0],
        meth.shape[1],
    )

    return meth


def prepare_methylation(meth, sites):
    """Subset, reorder, and fill absent GP-Age CpGs with NaN."""
    available = sites.intersection(meth.columns)
    missing = sites.difference(meth.columns)

    logger.info(
        "GP-Age CpGs: %d required, %d found, %d missing",
        len(sites),
        len(available),
        len(missing),
    )

    if len(available) == 0:
        raise ValueError(
            "None of the CpG sites required by GP-Age were found "
            "in the methylation dataset"
        )

    result = meth.reindex(columns=sites)

    return result


def impute_missing(meth, predictor, sites):
    """
    Impute missing methylation values using CpG-specific means
    from the GP-Age training data.
    """
    n_missing = int(meth.isna().sum().sum())

    if n_missing == 0:
        return meth

    logger.info("Imputing %d missing methylation values", n_missing)

    training = pd.DataFrame(
        np.asarray(predictor.X),
        columns=sites,
    )

    imputer = SimpleImputer(strategy="mean")
    imputer.fit(training)

    values = imputer.transform(meth)

    return pd.DataFrame(
        values,
        index=meth.index,
        columns=meth.columns,
    )


def predict_age(predictor, meth):
    """Predict biological age."""
    logger.info("Predicting biological age")

    pred_mean, pred_var = predictor.predict(meth.to_numpy(dtype=float))

    result = pd.DataFrame(
        {
            "GP_Age": pred_mean.ravel(),
            "GP_Age_SE": np.sqrt(pred_var.ravel()),
        },
        index=meth.index,
    )

    result.index.name = "Sample"

    return result


def read_age(infile, samples):
    """
    Read chronological ages.

    Accepts either:

        sample,age

    or a single-column file whose row order corresponds to the
    methylation samples.
    """
    age = pd.read_csv(infile, sep=None, engine="python")

    if age.shape[1] >= 2:
        sample_col = age.columns[0]
        age_col = age.columns[-1]

        age = age.set_index(sample_col)[age_col]
        age.index = age.index.astype(str)

        missing_samples = samples.difference(age.index)

        if len(missing_samples):
            raise ValueError(
                "Age information is missing for %d sample(s)"
                % len(missing_samples)
            )

        age = age.reindex(samples)

    else:
        age = age.iloc[:, 0]

        if len(age) != len(samples):
            raise ValueError(
                "Number of ages (%d) does not match number of samples (%d)"
                % (len(age), len(samples))
            )

        age.index = samples

    return pd.to_numeric(age, errors="raise")


def prediction_stats(predicted, observed):
    """Calculate age-prediction performance statistics."""
    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)

    rmse = np.sqrt(mean_squared_error(observed, predicted))
    medae = median_absolute_error(observed, predicted)
    mae = mean_absolute_error(observed, predicted)
    corr = np.corrcoef(observed, predicted)[0, 1]

    return pd.DataFrame(
        {
            "Statistic": ["RMSE", "MedAE", "MAE", "Pearson_r"],
            "Value": [rmse, medae, mae, corr],
        }
    )


def write_table(df, outfile=None, index=True):
    """Write a table to file or stdout."""
    if outfile:
        df.to_csv(outfile, sep="\t", index=index, float_format="%.4f")
        logger.info("Output saved to: %s", outfile)
    else:
        df.to_csv(sys.stdout, sep="\t", index=index, float_format="%.4f")


def run_gp_age(
    methylation_file,
    model="30",
    age_file=None,
    output_prefix=None,
    model_dir=None,
):
    """Run GP-Age prediction."""
    if model_dir is None:
        model_dir = get_model_dir()

    name = model_name(model)

    logger.info("GP-Age model: %s", name)

    sites = load_sites(model, model_dir)
    predictor = load_model(model, model_dir)

    meth = read_methylation(methylation_file)

    # Permit an "age" column embedded in the methylation matrix,
    # although a separate --age file is preferred.
    embedded_age = None

    if "age" in meth.columns:
        embedded_age = pd.to_numeric(meth.pop("age"), errors="coerce")

    meth = prepare_methylation(meth, sites)
    meth = impute_missing(meth, predictor, sites)

    predictions = predict_age(predictor, meth)

    if output_prefix:
        pred_file = "%s.GP_age.txt" % output_prefix
        write_table(predictions, pred_file)
    else:
        write_table(predictions)

    observed_age = None

    if age_file:
        observed_age = read_age(age_file, predictions.index)
    elif embedded_age is not None:
        observed_age = embedded_age.reindex(predictions.index)

    if observed_age is not None:
        valid = observed_age.notna()

        if valid.sum() < 2:
            logger.warning(
                "Fewer than two samples have valid chronological ages; "
                "performance statistics will not be calculated"
            )
            return predictions, None

        stats = prediction_stats(
            predictions.loc[valid, "GP_Age"],
            observed_age.loc[valid],
        )

        if output_prefix:
            stat_file = "%s.GP_age.stats.txt" % output_prefix
            write_table(stats, stat_file, index=False)
        else:
            sys.stderr.write("\nPrediction statistics:\n")
            stats.to_csv(
                sys.stderr,
                sep="\t",
                index=False,
                float_format="%.4f",
            )

        return predictions, stats

    return predictions, None


def main():
    parser = argparse.ArgumentParser(
        description="Predict biological age using the GP-Age model."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        metavar="FILE",
        help=(
            "DNA methylation matrix. Rows are CpG sites and columns "
            "are samples; the first column contains CpG IDs."
        ),
    )

    parser.add_argument(
        "-m",
        "--model",
        choices=MODEL_CHOICES,
        default="30",
        help="GP-Age model to use [default: %(default)s]",
    )

    parser.add_argument(
        "-a",
        "--age",
        metavar="FILE",
        help="Optional file containing chronological ages.",
    )

    parser.add_argument(
        "-o",
        "--output-prefix",
        metavar="PREFIX",
        help=(
            "Output prefix. Predictions are written to "
            "PREFIX.GP_age.txt and statistics to "
            "PREFIX.GP_age.stats.txt. If omitted, predictions "
            "are written to stdout."
        ),
    )

    parser.add_argument(
        "--model-dir",
        metavar="DIR",
        help="Directory containing GP-Age model and CpG files.",
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debugging information.",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    try:
        run_gp_age(
            methylation_file=args.input,
            model=args.model,
            age_file=args.age,
            output_prefix=args.output_prefix,
            model_dir=args.model_dir,
        )

    except (OSError, ValueError) as exc:
        logger.error("%s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

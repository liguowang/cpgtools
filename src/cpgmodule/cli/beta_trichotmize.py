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

"""Classify DNA methylation beta values into three methylation states.

For each sample, a three-component Bayesian Gaussian mixture model (BGMM) is
fitted to the beta-value distribution. Components are ordered by their means
and mapped to:

* 0: unmethylated
* 1: partially methylated
* 2: fully methylated
* -1: unassigned because the posterior probability is below the cutoff

The program writes one combined classification table and, optionally, a model
summary table.

Output files:

* PREFIX.methylation_states.tsv
    Each cell contains the assigned methylation state for one CpG in one sample.
    0 = Unmethylated
    1 = Partially methylated
    2 = Fully methylated
    -1 = Unassigned (maximum posterior probability below --prob_cutoff
                     
* PREFIX.assignment_probability.ts
  The posterior probability associated with the assigned methylation state.
  Values range from 0 to 1.
  Higher values indicate greater confidence in the state assignment.
  
* PREFIX.probability_0_unmethylated.tsv
  Posterior probability that each CpG belongs to the unmethylated component.
  
* PREFIX.probability_1_partially_methylated.tsv
  Posterior probability that each CpG belongs to the partially methylated component.
  
* PREFIX.probability_2_fully_methylated.tsv
  Posterior probability that each CpG belongs to the fully methylated component.

* PREFIX.bgmm_summary.tsv (Optional)
  Contains one row per sample summarizing the fitted Bayesian Gaussian mixture model.

"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.mixture import BayesianGaussianMixture

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


STATE_LABELS = {
    0: "unmethylated",
    1: "partially_methylated",
    2: "fully_methylated",
}


class TrichotomizeError(RuntimeError):
    """Raised when BGMM fitting or classification fails."""


def build_parser() -> argparse.ArgumentParser:
    """Build and return the command-line parser."""
    parser = argparse.ArgumentParser(
        prog="cpgtools beta_trichotomize",
        description=(
            "Classify beta values into unmethylated, partially methylated, "
            "and fully methylated states using per-sample Bayesian Gaussian "
            "mixture models."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        type=Path,
        help=(
            "Input beta-value matrix. The first column contains CpG IDs and "
            "remaining columns contain samples. Compressed input is supported."
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
        "-c",
        "--prob_cutoff",
        "--prob-cut",
        dest="prob_cutoff",
        type=float,
        default=0.95,
        help="Minimum posterior probability required for state assignment.",
    )
    parser.add_argument(
        "-s",
        "--seed",
        dest="random_state",
        type=int,
        default=99,
        help="Random seed used for BGMM fitting.",
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=5000,
        help="Maximum number of BGMM iterations per sample.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-3,
        help="BGMM convergence tolerance.",
    )
    parser.add_argument(
        "--weight_concentration_prior",
        type=float,
        default=None,
        help=(
            "Optional Dirichlet-process weight concentration prior. "
            "By default, scikit-learn chooses the prior."
        ),
    )
    parser.add_argument(
        "--covariance_type",
        choices=("full", "tied", "diag", "spherical"),
        default="full",
        help="BGMM covariance parameterization.",
    )
    parser.add_argument(
        "--na_policy",
        choices=("drop", "per_sample"),
        default="per_sample",
        help=(
            "'drop' removes CpGs with any missing value before fitting all "
            "samples; 'per_sample' fits each sample using its available values."
        ),
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Write a model-summary table.",
    )
    parser.add_argument(
        "--long_format",
        action="store_true",
        help=(
            "Write classification output in long format. By default, separate "
            "state and probability matrices are written."
        ),
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
        parser.error(f"Input file does not exist: {args.input_file}")
    if not 0.0 <= args.prob_cutoff <= 1.0:
        parser.error("--prob_cutoff must be between 0 and 1")
    if args.max_iter < 1:
        parser.error("--max_iter must be at least 1")
    if args.tol <= 0:
        parser.error("--tol must be greater than zero")
    if (
        args.weight_concentration_prior is not None
        and args.weight_concentration_prior <= 0
    ):
        parser.error("--weight_concentration_prior must be greater than zero")


def read_beta_matrix(path: Path, na_policy: str) -> pd.DataFrame:
    """Read and validate the input beta-value matrix."""
    printlog(f'Reading beta-value matrix: "{path}"')

    try:
        data = pd.read_csv(
            path,
            sep=None,
            engine="python",
            compression="infer",
            index_col=0,
        )
    except Exception as exc:
        raise TrichotomizeError(f"Cannot read input file {path}: {exc}") from exc

    if data.shape[1] == 0:
        raise TrichotomizeError("Input contains no sample columns.")

    data.index = data.index.astype(str)
    data.columns = data.columns.astype(str)

    if data.index.duplicated().any():
        duplicates = data.index[data.index.duplicated()].unique().tolist()
        raise TrichotomizeError(
            "Duplicate CpG IDs were found: " + ", ".join(duplicates[:10])
        )

    if data.columns.duplicated().any():
        duplicates = data.columns[data.columns.duplicated()].tolist()
        raise TrichotomizeError(
            "Duplicate sample IDs were found: " + ", ".join(duplicates)
        )

    data = data.apply(pd.to_numeric, errors="coerce")

    finite = data.to_numpy(dtype=float)
    out_of_range_mask = np.isfinite(finite) & ((finite < 0.0) | (finite > 1.0))
    if out_of_range_mask.any():
        raise TrichotomizeError(
            f"Input contains {int(out_of_range_mask.sum())} finite beta values "
            "outside [0, 1]."
        )

    if na_policy == "drop":
        before = len(data)
        data = data.dropna(axis=0, how="any")
        printlog(f"Removed {before - len(data)} CpGs containing missing values.")

    if data.empty:
        raise TrichotomizeError("No CpGs remain after filtering.")

    printlog(f"Total samples: {data.shape[1]}")
    printlog(f"Total CpGs: {data.shape[0]}")
    return data


def fit_sample(
    sample_name: str,
    values: pd.Series,
    *,
    random_state: int,
    max_iter: int,
    tol: float,
    covariance_type: str,
    weight_concentration_prior: Optional[float],
) -> tuple[BayesianGaussianMixture, np.ndarray, np.ndarray]:
    """Fit a three-component BGMM and return ordered state probabilities."""
    valid_mask = values.notna().to_numpy()
    x = values.to_numpy(dtype=float)[valid_mask]

    if x.size < 3:
        raise TrichotomizeError(
            f'Sample "{sample_name}" has fewer than three valid beta values.'
        )

    unique_values = np.unique(x)
    if unique_values.size < 3:
        raise TrichotomizeError(
            f'Sample "{sample_name}" has fewer than three distinct beta values.'
        )

    kwargs = {
        "n_components": 3,
        "covariance_type": covariance_type,
        "max_iter": max_iter,
        "tol": tol,
        "random_state": random_state,
        "weight_concentration_prior_type": "dirichlet_process",
        "n_init": 1,
    }
    if weight_concentration_prior is not None:
        kwargs["weight_concentration_prior"] = weight_concentration_prior

    model = BayesianGaussianMixture(**kwargs)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model.fit(x.reshape(-1, 1))

    for warning in caught:
        printlog(f'Warning for sample "{sample_name}": {warning.message}')

    means = model.means_.reshape(-1)
    state_to_component = np.argsort(means)
    component_to_state = np.empty(3, dtype=int)
    component_to_state[state_to_component] = np.arange(3)

    probabilities = np.full((len(values), 3), np.nan, dtype=float)
    raw_probabilities = model.predict_proba(x.reshape(-1, 1))
    probabilities[valid_mask] = raw_probabilities[:, state_to_component]

    return model, probabilities, component_to_state


def classify_probabilities(
    probabilities: np.ndarray,
    probability_cutoff: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ordered state probabilities into labels and confidence values."""
    states = np.full(probabilities.shape[0], -1, dtype=int)
    confidence = np.full(probabilities.shape[0], np.nan, dtype=float)

    valid = np.all(np.isfinite(probabilities), axis=1)
    if not np.any(valid):
        return states, confidence

    best_state = np.argmax(probabilities[valid], axis=1)
    best_probability = np.max(probabilities[valid], axis=1)

    assigned = best_probability >= probability_cutoff
    states_valid = np.where(assigned, best_state, -1)

    states[valid] = states_valid
    confidence[valid] = best_probability
    return states, confidence


def process_samples(
    data: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    dict[str, np.ndarray],
    list[dict[str, object]],
]:
    """Fit all samples and collect classifications and model summaries."""
    state_results: dict[str, np.ndarray] = {}
    confidence_results: dict[str, np.ndarray] = {}
    probability_results: dict[str, np.ndarray] = {}
    summaries: list[dict[str, object]] = []

    for sample_name in data.columns:
        printlog(f'Fitting BGMM for sample: "{sample_name}"')

        model, ordered_probabilities, component_to_state = fit_sample(
            sample_name,
            data[sample_name],
            random_state=args.random_state,
            max_iter=args.max_iter,
            tol=args.tol,
            covariance_type=args.covariance_type,
            weight_concentration_prior=args.weight_concentration_prior,
        )

        states, confidence = classify_probabilities(
            ordered_probabilities,
            args.prob_cutoff,
        )

        state_results[sample_name] = states
        confidence_results[sample_name] = confidence
        probability_results[sample_name] = ordered_probabilities

        component_order = np.argsort(model.means_.reshape(-1))
        ordered_means = model.means_.reshape(-1)[component_order]
        ordered_weights = model.weights_.reshape(-1)[component_order]

        summaries.append(
            {
                "sample": sample_name,
                "unmethylated_mean": ordered_means[0],
                "partially_methylated_mean": ordered_means[1],
                "fully_methylated_mean": ordered_means[2],
                "unmethylated_weight": ordered_weights[0],
                "partially_methylated_weight": ordered_weights[1],
                "fully_methylated_weight": ordered_weights[2],
                "converged": bool(model.converged_),
                "n_iter": int(model.n_iter_),
                "lower_bound": float(model.lower_bound_),
                "valid_cpgs": int(data[sample_name].notna().sum()),
                "assigned_cpgs": int(np.sum(states >= 0)),
                "unassigned_cpgs": int(np.sum(states < 0)),
            }
        )

    return (
        state_results,
        confidence_results,
        probability_results,
        summaries,
    )


def write_matrix_outputs(
    data: pd.DataFrame,
    state_results: dict[str, np.ndarray],
    confidence_results: dict[str, np.ndarray],
    probability_results: dict[str, np.ndarray],
    out_prefix: Path,
) -> None:
    """Write state, confidence, and class-probability matrices."""
    states = pd.DataFrame(state_results, index=data.index)
    states.index.name = data.index.name or "CpG_ID"

    confidence = pd.DataFrame(confidence_results, index=data.index)
    confidence.index.name = states.index.name

    state_path = Path(f"{out_prefix}.methylation_states.tsv")
    confidence_path = Path(f"{out_prefix}.assignment_probability.tsv")
    printlog(f'Writing methylation-state matrix: "{state_path}"')
    states.to_csv(state_path, sep="\t", na_rep="-1")
    printlog(f'Writing assignment-probability matrix: "{confidence_path}"')
    confidence.to_csv(
        confidence_path,
        sep="\t",
        na_rep="NA",
        float_format="%.6f",
    )

    for state, state_name in STATE_LABELS.items():
        probability_table = pd.DataFrame(
            {
                sample: probabilities[:, state]
                for sample, probabilities in probability_results.items()
            },
            index=data.index,
        )
        probability_table.index.name = states.index.name
        probability_path = Path(
            f"{out_prefix}.probability_{state}_{state_name}.tsv"
        )
        printlog(f'Writing state-probability matrix: "{probability_path}"')
        probability_table.to_csv(
            probability_path,
            sep="\t",
            na_rep="NA",
            float_format="%.6f",
        )


def write_long_output(
    data: pd.DataFrame,
    state_results: dict[str, np.ndarray],
    confidence_results: dict[str, np.ndarray],
    probability_results: dict[str, np.ndarray],
    out_prefix: Path,
) -> None:
    """Write one row per CpG and sample."""
    records = []

    for sample in data.columns:
        probabilities = probability_results[sample]
        states = state_results[sample]
        confidence = confidence_results[sample]

        sample_frame = pd.DataFrame(
            {
                "CpG_ID": data.index,
                "sample": sample,
                "beta_value": data[sample].to_numpy(dtype=float),
                "probability_0": probabilities[:, 0],
                "probability_1": probabilities[:, 1],
                "probability_2": probabilities[:, 2],
                "assigned_state": states,
                "assignment_probability": confidence,
            }
        )
        records.append(sample_frame)

    result = pd.concat(records, ignore_index=True)
    output_path = Path(f"{out_prefix}.trichotomized.long.tsv")
    printlog(f'Writing long-format classifications: "{output_path}"')
    result.to_csv(
        output_path,
        sep="\t",
        index=False,
        na_rep="NA",
        float_format="%.6f",
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args, parser)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        data = read_beta_matrix(args.input_file, args.na_policy)

        (
            state_results,
            confidence_results,
            probability_results,
            summaries,
        ) = process_samples(data, args)

        if args.long_format:
            write_long_output(
                data,
                state_results,
                confidence_results,
                probability_results,
                args.out_prefix,
            )
        else:
            write_matrix_outputs(
                data,
                state_results,
                confidence_results,
                probability_results,
                args.out_prefix,
            )

        if args.report:
            summary_path = Path(f"{args.out_prefix}.bgmm_summary.tsv")
            printlog(f'Writing BGMM summary: "{summary_path}"')
            pd.DataFrame(summaries).to_csv(
                summary_path,
                sep="\t",
                index=False,
                float_format="%.6f",
            )

    except TrichotomizeError as exc:
        parser.exit(1, f"Error: {exc}\n")
    except Exception as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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
Estimate how different the means of two groups are using a Bayesian approach.

This command uses a Metropolis-Hastings Markov chain Monte Carlo (MCMC)
procedure inspired by John Kruschke's BEST framework (Bayesian Estimation
Supersedes the t test). For each CpG/probe, it estimates:

* posterior mean for group 1
* posterior mean for group 2
* posterior difference of means
* 95% highest density interval (HDI) for the difference
* posterior probability that the difference is on the same side of zero
  as the posterior median

Notes
-----
This method is substantially slower than a standard t test because an MCMC
chain is fitted independently for each CpG/probe. Multiple worker processes
are therefore recommended for large datasets.
"""

import argparse
import collections
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import stats

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog, read_grp_file1


def dt(x, mu, sig):
    """Return Student-t log-density values for observations in ``x``."""
    x = np.asarray(x, dtype=float)
    if sig <= 0 or x.size < 2:
        return np.full(x.shape, -np.inf, dtype=float)
    return stats.t.logpdf(x, loc=mu, scale=sig, df=x.size - 1)


def dnorm(x, mu, sig):
    """Return the normal log-density."""
    if sig <= 0:
        return -np.inf
    return stats.norm.logpdf(x, loc=mu, scale=sig)


def dexp(x, rate):
    """Return the exponential log-density using a rate parameter."""
    if x <= 0 or rate <= 0:
        return -np.inf
    return stats.expon.logpdf(x, scale=1.0 / rate)


def likelihood(s1, s2, params):
    """Return the log-likelihood for two groups."""
    mu1, sig1, mu2, sig2 = params
    return np.sum(dt(s1, mu1, sig1)) + np.sum(dt(s2, mu2, sig2))


def prior(s1, s2, params):
    """Return the log-prior density for model parameters."""
    mu1, sig1, mu2, sig2 = params
    pooled = np.concatenate((s1, s2))

    prior_mean = float(np.mean(pooled))
    pooled_std = float(np.std(pooled, ddof=0))

    # Preserve the historical very broad prior while avoiding a zero-width
    # normal prior when all pooled values are identical.
    prior_std = max(1000.0 * pooled_std, np.finfo(float).eps)

    return (
        dnorm(mu1, prior_mean, prior_std)
        + dnorm(mu2, prior_mean, prior_std)
        + dexp(sig1, 0.1)
        + dexp(sig2, 0.1)
    )


def posterior(s1, s2, params):
    """Return the log posterior density."""
    return likelihood(s1, s2, params) + prior(s1, s2, params)


def compute_hdi(chain, interval=0.95):
    """
    Compute the highest density interval (HDI).

    Parameters
    ----------
    chain : array-like
        Posterior samples.
    interval : float, optional
        Credible mass. Default is 0.95.

    Returns
    -------
    tuple[float, float]
        Lower and upper HDI bounds.
    """
    samples = np.asarray(chain, dtype=float)

    if samples.size == 0:
        raise ValueError("Cannot compute HDI from an empty chain")

    if not 0 < interval < 1:
        raise ValueError("HDI interval must be between 0 and 1")

    samples = np.sort(samples)

    n_samples = samples.size
    n_credible = int(np.floor(interval * n_samples))

    if n_credible < 1:
        raise ValueError("Not enough posterior samples to compute HDI")

    if n_credible >= n_samples:
        return float(samples[0]), float(samples[-1])

    widths = samples[n_credible:] - samples[: n_samples - n_credible]
    best = int(np.argmin(widths))

    return (
        float(samples[best]),
        float(samples[best + n_credible]),
    )


def _initial_scale(values):
    """Return a positive initial scale for one group."""
    if values.size < 2:
        return np.finfo(float).eps

    scale = float(np.std(values, ddof=0))
    return max(scale, np.finfo(float).eps)


def beta_bayes(
    probe_id,
    s1,
    s2,
    seed,
    niter=5000,
    nburn_in=500,
):
    """
    Fit the Bayesian two-group model for one CpG/probe.

    Returns
    -------
    tuple
        probe ID, posterior means, difference, HDI bounds, probability.
    """
    s1 = np.asarray(s1, dtype=float)
    s2 = np.asarray(s2, dtype=float)

    if s1.size < 2 or s2.size < 2:
        return (
            probe_id,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )

    if niter <= 1:
        raise ValueError("niter must be greater than 1")

    if nburn_in < 0 or nburn_in >= niter:
        raise ValueError("burn-in must be >= 0 and smaller than niter")

    rng = np.random.default_rng(seed)

    sig1 = _initial_scale(s1)
    sig2 = _initial_scale(s2)

    parameters = np.array(
        [
            float(np.mean(s1)),
            sig1,
            float(np.mean(s2)),
            sig2,
        ],
        dtype=float,
    )

    increment = max(
        (sig1 + sig2) / 10.0,
        np.finfo(float).eps,
    )

    mu1_samples = []
    mu2_samples = []

    current_log_post = posterior(s1, s2, parameters)

    for iteration in range(1, niter):
        candidate = parameters + rng.normal(0.0, increment, 4)

        if candidate[1] <= 0 or candidate[3] <= 0:
            if iteration >= nburn_in:
                mu1_samples.append(parameters[0])
                mu2_samples.append(parameters[2])
            continue

        candidate_log_post = posterior(s1, s2, candidate)
        log_alpha = candidate_log_post - current_log_post

        if math.log(rng.uniform()) < min(0.0, log_alpha):
            parameters = candidate
            current_log_post = candidate_log_post

        if iteration >= nburn_in:
            mu1_samples.append(parameters[0])
            mu2_samples.append(parameters[2])

    mu1_samples = np.asarray(mu1_samples, dtype=float)
    mu2_samples = np.asarray(mu2_samples, dtype=float)

    if mu1_samples.size == 0:
        return (
            probe_id,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )

    est_mu1 = float(np.mean(mu1_samples))
    est_mu2 = float(np.mean(mu2_samples))

    diff = mu1_samples - mu2_samples
    diff_median = float(np.median(diff))

    if diff_median < 0:
        probability = float(np.mean(diff < 0))
    elif diff_median > 0:
        probability = float(np.mean(diff > 0))
    else:
        probability = 0.5

    hdi_low, hdi_high = compute_hdi(diff, interval=0.95)

    return (
        probe_id,
        est_mu1,
        est_mu2,
        est_mu1 - est_mu2,
        hdi_low,
        hdi_high,
        probability,
    )


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
            "Data matrix containing sample IDs in the first row and CpG/probe "
            "IDs in the first column. Non-numeric values in the data matrix "
            "are treated as missing and ignored. Plain text, .gz, and .bz2 "
            "inputs are supported."
        ),
    )

    parser.add_argument(
        "-g",
        "--group",
        dest="group_file",
        required=True,
        help=(
            "Comma-separated two-column group file with a header. Column 1 "
            "contains sample IDs and column 2 contains group IDs. Exactly two "
            "groups are required."
        ),
    )

    parser.add_argument(
        "-n",
        "--niter",
        dest="n_iter",
        type=int,
        default=5000,
        help="Number of MCMC iterations [default: %(default)s].",
    )

    parser.add_argument(
        "-b",
        "--burnin",
        dest="n_burn",
        type=int,
        default=500,
        help="Number of initial MCMC samples to discard [default: %(default)s].",
    )

    parser.add_argument(
        "-p",
        "--processor",
        dest="n_process",
        type=int,
        default=1,
        help="Number of worker processes [default: %(default)s].",
    )

    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=99,
        help="Random-number seed [default: %(default)s].",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="out_file",
        required=True,
        help="Output prefix. Results are written to <prefix>.bayes.tsv.",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def validate_args(args, parser=None):
    """Validate parsed arguments."""
    errors = []

    if args.n_iter <= 1:
        errors.append("--niter must be greater than 1")

    if args.n_burn < 0:
        errors.append("--burnin must be zero or greater")

    if args.n_burn >= args.n_iter:
        errors.append("--burnin must be smaller than --niter")

    if args.n_process <= 0:
        errors.append("--processor must be a positive integer")

    if errors:
        message = "; ".join(errors)
        if parser is not None:
            parser.error(message)
        raise ValueError(message)


def read_groups(group_file):
    """
    Read the group file and validate that exactly two groups are present.

    Returns
    -------
    tuple
        sample-to-group mapping, ordered group IDs, group-to-samples mapping.
    """
    printlog(f'Read group file "{group_file}" ...')
    samples, groups = read_grp_file1(group_file)

    sample_to_group = dict(zip(samples, groups))

    group_to_samples = collections.defaultdict(list)
    for sample, group in zip(samples, groups):
        group_to_samples[group].append(sample)

    group_ids = sorted(group_to_samples)

    for group in group_ids:
        print(
            f"\tGroup {group} has {len(group_to_samples[group])} samples:",
            file=sys.stderr,
        )
        print(
            "\t\t" + ",".join(group_to_samples[group]),
            file=sys.stderr,
        )

    if len(group_ids) != 2:
        raise ValueError("Exactly two groups are required")

    return sample_to_group, group_ids, group_to_samples


def parse_matrix(input_file, sample_to_group, group_ids):
    """
    Yield Bayesian jobs from an input beta-value matrix.

    Each yielded item contains:
        probe_id, group1_values, group2_values
    """
    printlog(f'Read data file "{input_file}" ...')

    sample_group_ids = None

    for line_number, raw_line in enumerate(
        ireader.reader(input_file),
        start=1,
    ):
        fields = raw_line.split()

        if not fields:
            continue

        if line_number == 1:
            sample_ids = fields[1:]

            missing = [
                sample
                for sample in sample_to_group
                if sample not in sample_ids
            ]
            if missing:
                raise ValueError(
                    "Cannot find sample ID(s) in data file: "
                    + ", ".join(missing)
                )

            sample_group_ids = [
                sample_to_group.get(sample)
                for sample in sample_ids
            ]
            continue

        if sample_group_ids is None:
            raise ValueError("Input matrix does not contain a valid header")

        probe_id = fields[0]
        beta_values = fields[1:]

        group1 = []
        group2 = []

        for group_id, value in zip(sample_group_ids, beta_values):
            try:
                beta = float(value)
            except ValueError:
                continue

            if not np.isfinite(beta):
                continue

            if group_id == group_ids[0]:
                group1.append(beta)
            elif group_id == group_ids[1]:
                group2.append(beta)

        yield (
            probe_id,
            np.asarray(group1, dtype=float),
            np.asarray(group2, dtype=float),
        )


def _worker(args):
    """Multiprocessing wrapper for one Bayesian model fit."""
    return beta_bayes(*args)


def run_bayes(
    input_file,
    group_file,
    output_prefix,
    n_iter=5000,
    n_burn=500,
    n_process=1,
    seed=99,
):
    """
    Run Bayesian differential methylation analysis.

    Returns
    -------
    str
        Path to the generated TSV file.
    """
    sample_to_group, group_ids, _ = read_groups(group_file)

    jobs = []
    for index, (probe_id, group1, group2) in enumerate(
        parse_matrix(
            input_file,
            sample_to_group,
            group_ids,
        )
    ):
        # Derive a deterministic per-probe seed so independent workers do not
        # all start from the same random stream.
        job_seed = int(seed + index)

        jobs.append(
            (
                probe_id,
                group1,
                group2,
                job_seed,
                n_iter,
                n_burn,
            )
        )

    output_path = Path(f"{output_prefix}.bayes.tsv")
    if output_path.parent != Path("."):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    if n_process == 1:
        results = [_worker(job) for job in jobs]
    else:
        with Pool(processes=n_process) as pool:
            results = pool.map(_worker, jobs)

    with open(output_path, "w") as fout:
        print(
            "\t".join(
                [
                    "ID",
                    "mu1",
                    "mu2",
                    "mu_diff",
                    "mu_diff (95% HDI)",
                    "Probability",
                ]
            ),
            file=fout,
        )

        for result in results:
            (
                probe_id,
                mu1,
                mu2,
                mu_diff,
                hdi_low,
                hdi_high,
                probability,
            ) = result

            print(
                f"{probe_id}\t"
                f"{mu1:.6f}\t"
                f"{mu2:.6f}\t"
                f"{mu_diff:.6f}\t"
                f"({hdi_low:.6f},{hdi_high:.6f})\t"
                f"{probability:.6f}",
                file=fout,
            )

    return str(output_path)


def main(argv=None):
    """
    Command-line entry point.

    Parameters
    ----------
    argv : list[str] or None
        Optional command-line argument list. When None, argparse reads
        ``sys.argv``. Passing a list makes this function directly callable
        from the CpGtools dispatcher and from tests.

    Returns
    -------
    int
        Process-style return code.
    """
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        validate_args(args, parser=parser)

        run_bayes(
            input_file=args.input_file,
            group_file=args.group_file,
            output_prefix=args.out_file,
            n_iter=args.n_iter,
            n_burn=args.n_burn,
            n_process=args.n_process,
            seed=args.seed,
        )
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())

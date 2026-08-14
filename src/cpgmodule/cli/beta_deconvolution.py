#!/usr/bin/env python3
"""Estimate cell/tissue fractions from a DNA methylation beta-value matrix.

The input matrix must contain CpG identifiers in the first column and one or
more sample beta-value columns. Reference atlases are bundled with CpGtools
as package data.
"""

from __future__ import annotations

import argparse
import math
import multiprocessing as mp
import sys
from pathlib import Path
from typing import Optional, Sequence, Tuple

import matplotlib
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, nnls
from importlib.resources import files


ATLAS_DIR = files("cpgmodule") / "deconv_data"

ATLAS_FILES = {
    "CellTrace": "Reference_CellTrace.beta.tsv",
    "EpiDISH": "Reference_EpiDISH.beta.tsv",
}

NR_SAMPLE_XTICKS = 30
FIG_SIZE = (15, 7)
COLOR_MAP = "tab10"
OTHERS_THRESHOLD = 0.01

# Optional collapse of immune-cell subtypes. It is applied only when the
# corresponding labels are present in the selected reference atlas.
IMMUNE_CELL_MAPPING = {
    "Bcell.pred": ("Bmem.pred", "Bnv.pred"),
    "CD4T.pred": ("CD4mem.pred", "CD4nv.pred"),
    "CD8T.pred": ("CD8mem.pred", "CD8nv.pred"),
    "Granulocyte.pred": ("Neu.pred", "Eos.pred", "Bas.pred"),
}


class DeconvolutionError(RuntimeError):
    """Raised when a deconvolution input or fit is invalid."""


def read_beta_matrix(path: Path) -> pd.DataFrame:
    """Read a tabular beta-value matrix and standardize its CpG column."""
    if not path.is_file():
        raise DeconvolutionError(f"Input file does not exist: {path}")

    try:
        data = pd.read_csv(path, sep=None, engine="python", compression="infer")
    except Exception as exc:  # pragma: no cover - pandas supplies details
        raise DeconvolutionError(f"Cannot read {path}: {exc}") from exc

    if data.shape[1] < 2:
        raise DeconvolutionError(
            f"{path} must contain a CpG column and at least one value column."
        )

    data = data.rename(columns={data.columns[0]: "acc"})
    data["acc"] = data["acc"].astype(str)
    data = data.drop_duplicates(subset="acc", keep="first")
    return data


def validate_numeric_columns(data: pd.DataFrame, label: str) -> pd.DataFrame:
    """Convert all non-CpG columns to numeric values."""
    result = data.copy()
    for column in result.columns[1:]:
        result[column] = pd.to_numeric(result[column], errors="coerce")

    if result.iloc[:, 1:].notna().sum().sum() == 0:
        raise DeconvolutionError(f"No numeric beta values were found in {label}.")
    return result


def hide_small_components(data: pd.DataFrame) -> pd.DataFrame:
    """Collapse fractions below the plotting threshold into an 'other' row."""
    small = data.where(data < OTHERS_THRESHOLD, 0.0).sum(axis=0)
    cleaned = data.mask(data < OTHERS_THRESHOLD, 0.0)
    return pd.concat([cleaned, pd.DataFrame([small], index=["other"])])


def generate_bar_styles(number_components: int):
    """Return color/hatch pairs for a stacked bar plot."""
    matplotlib.rcParams["hatch.linewidth"] = 0.3
    hatch_types = [None, "xxx", "...", "O", "++"]
    hatches = hatch_types[: max(1, math.ceil(number_components / 7))]
    number_colors = max(1, math.ceil(number_components / len(hatches)))

    cmap = matplotlib.colormaps.get_cmap(COLOR_MAP)
    norm = matplotlib.colors.Normalize(vmin=0.0, vmax=float(number_colors))
    colors = [cmap(norm(index)) for index in range(number_colors)]

    styles = []
    for index in range(number_components):
        color_index = index % number_colors
        hatch_index = min(index // number_colors, len(hatches) - 1)
        styles.append((colors[color_index], hatches[hatch_index]))
    return styles


def plot_results(data: pd.DataFrame, output_path: Path, show: bool = False) -> None:
    """Create a stacked bar plot of estimated fractions."""
    plot_data = hide_small_components(data)
    number_components, number_samples = plot_data.shape
    styles = generate_bar_styles(number_components)

    figure, axis = plt.subplots(figsize=FIG_SIZE)
    sample_positions = np.arange(number_samples)
    bottom = np.zeros(number_samples, dtype=float)

    for index in range(number_components):
        values = plot_data.iloc[index].to_numpy(dtype=float)
        axis.bar(
            sample_positions,
            values,
            edgecolor="white",
            width=0.85,
            label=plot_data.index[index],
            bottom=bottom,
            color=styles[index][0],
            hatch=styles[index][1],
        )
        bottom += values

    axis.set_xticks(sample_positions)
    axis.set_xticklabels(
        [name[:NR_SAMPLE_XTICKS] for name in plot_data.columns],
        rotation="vertical",
        fontsize=9,
    )
    axis.set_xlabel("Sample")
    axis.set_ylabel("Estimated fraction")
    axis.set_ylim(0.0, 1.0)
    axis.set_xlim(-0.6, number_samples - 0.4)
    axis.legend(loc="upper left", bbox_to_anchor=(1, 1), ncol=1)
    axis.set_title(f"DNA methylation deconvolution\n{output_path.name}")
    figure.tight_layout(rect=[0, 0, 0.83, 1])
    figure.savefig(output_path, dpi=300)

    if show:
        plt.show()
    plt.close(figure)


def fit_sample(job):
    """Fit one sample and return fractions, L2 residual, RMSE, and fit status.

    Parameters
    ----------
    job : tuple
        ``(sample_name, reference_matrix, sample_values, method)``.

    Returns
    -------
    tuple
        ``(sample_name, mixture, l2, rmse, n_cpg, error)``. On failure,
        ``mixture``, ``l2``, and ``rmse`` are ``None`` and ``error`` contains
        the failure message.
    """
    sample_name, reference_matrix, sample_values, method = job

    try:
        X = np.asarray(reference_matrix, dtype=float)
        y = np.asarray(sample_values, dtype=float)

        valid = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        X = X[valid]
        y = y[valid]
        n_cpg = int(y.size)

        if n_cpg == 0:
            raise DeconvolutionError("no complete CpGs are available for fitting")
        if X.shape[1] == 0:
            raise DeconvolutionError("the reference atlas has no cell-type columns")

        if method == "nnls":
            mixture, _ = nnls(X, y, maxiter=200)
        elif method == "bounded":
            result = lsq_linear(X, y, bounds=(0.0, 1.0), method="trf")
            if not result.success:
                raise DeconvolutionError(result.message)
            mixture = result.x
        else:  # Defensive check; argparse already restricts this value.
            raise DeconvolutionError(f"unsupported fitting method: {method}")

        mixture = np.clip(mixture, 0.0, None)
        total = float(mixture.sum())
        if total <= 0.0:
            raise DeconvolutionError("estimated fractions sum to zero")
        mixture /= total

        predicted = X @ mixture
        residual_vector = y - predicted
        l2 = float(np.linalg.norm(residual_vector, ord=2))
        rmse = float(np.sqrt(np.mean(np.square(residual_vector))))

        return sample_name, mixture, l2, rmse, n_cpg, None
    except Exception as exc:
        return sample_name, None, None, None, 0, str(exc)


class Deconvolve:
    """Run reference-based deconvolution of methylation beta values."""

    def __init__(
        self,
        atlas_path: Path,
        samples_path: Path,
        output_prefix: Optional[Path],
        method: str = "nnls",
        processes: Optional[int] = None,
        residuals: bool = False,
        slim: bool = False,
        plot: bool = False,
    ) -> None:
        self.atlas_path = atlas_path
        self.samples_path = samples_path
        
        self.method = method
        self.processes = processes
        self.write_residuals = residuals
        self.slim = slim
        self.make_plot = plot
        
        if output_prefix is None:
            self.output_prefix = Path(self._sample_stem(samples_path))
        else:
            self.output_prefix = output_prefix
        
        self.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        self.atlas, self.samples = self._load_inputs()

    @staticmethod
    def _sample_stem(path: Path) -> str:
        name = path.name
        if name.endswith(".gz"):
            name = name[:-3]
        return Path(name).stem

    def _load_inputs(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        atlas = validate_numeric_columns(read_beta_matrix(self.atlas_path), "atlas")
        samples = validate_numeric_columns(
            read_beta_matrix(self.samples_path), "sample matrix"
        )

        shared_cpgs = atlas[["acc"]].merge(samples[["acc"]], on="acc", how="inner")
        if shared_cpgs.empty:
            raise DeconvolutionError(
                "The sample matrix and reference atlas have no CpGs in common."
            )

        atlas = shared_cpgs.merge(atlas, on="acc", how="left")
        samples = shared_cpgs.merge(samples, on="acc", how="left")
        return atlas, samples

    def run(self) -> None:
        component_names = list(self.atlas.columns[1:])
        sample_names = list(self.samples.columns[1:])
        reference = self.atlas.iloc[:, 1:].to_numpy(dtype=float)

        jobs = [
            (
                sample_name,
                reference,
                self.samples[sample_name].to_numpy(dtype=float),
                self.method,
            )
            for sample_name in sample_names
        ]

        if self.processes == 1 or len(jobs) == 1:
            fitted = [fit_sample(job) for job in jobs]
        else:
            with mp.Pool(processes=self.processes) as pool:
                fitted = pool.map(fit_sample, jobs)

        mixtures = {}
        fit_metrics = {}
        failures = []
        for sample_name, mixture, l2, rmse, n_cpg, error in fitted:
            if error is not None:
                failures.append(f"{sample_name}: {error}")
                continue
            mixtures[sample_name] = mixture
            fit_metrics[sample_name] = {
                "residual_l2": l2,
                "rmse": rmse,
                "n_cpg": n_cpg,
            }

        if not mixtures:
            details = "; ".join(failures) if failures else "unknown fitting error"
            raise DeconvolutionError(f"No samples could be deconvolved: {details}")

        if failures:
            print(
                "Warning: some samples were skipped: " + "; ".join(failures),
                file=sys.stderr,
            )

        full_results = pd.DataFrame(mixtures, index=component_names)
        full_path = Path(f"{self.output_prefix}_deconv.tsv")
        self._write_table(full_results.T, full_path)

        plot_data = full_results
        combined_results = self._combine_immune_subtypes(full_results)
        if not combined_results.equals(full_results):
            combined_path = Path(f"{self.output_prefix}_combined_deconv.tsv")
            self._write_table(combined_results.T, combined_path)
            plot_data = combined_results

        if self.write_residuals:
            metrics_path = Path(f"{self.output_prefix}_fit_metrics.tsv")
            metrics_data = pd.DataFrame.from_dict(fit_metrics, orient="index")
            metrics_data.index.name = "sample"
            metrics_data.to_csv(
                metrics_path,
                sep="\t",
                index=True,
                float_format="%.6g",
            )

        if self.make_plot:
            plot_path = Path(f"{self.output_prefix}_deconv_plot.png")
            plot_results(plot_data, plot_path)

    def _write_table(self, data: pd.DataFrame, path: Path) -> None:
        data.index.name = "sample"
        data.to_csv(
            path,
            sep="\t",
            index=not self.slim,
            header=not self.slim,
            float_format="%.6f",
        )

    @staticmethod
    def _combine_immune_subtypes(data: pd.DataFrame) -> pd.DataFrame:
        combined = data.copy()
        for broad_type, subtypes in IMMUNE_CELL_MAPPING.items():
            existing = [subtype for subtype in subtypes if subtype in combined.index]
            if not existing:
                continue
            combined.loc[broad_type] = combined.loc[existing].sum(axis=0)
            combined = combined.drop(index=existing)
        return combined


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate cell or tissue fractions from a DNA methylation beta-value "
            "matrix using a bundled reference atlas."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "samples_path",
        type=Path,
        help="Input beta-value matrix; CpGs in column 1 and samples in later columns.",
    )
    parser.add_argument(
        "--atlas",
        choices=ATLAS_FILES.keys(),
        default="CellTrace",
        help="Reference methylation atlas.",
    )
    parser.add_argument(
        "--method",
        choices=("nnls", "bounded"),
        default="nnls",
        help=(
            "Fitting method. 'nnls' estimates non-negative coefficients and "
            "normalizes them to sum to one after fitting. 'bounded' constrains "
            "coefficients to the interval [0, 1] during fitting and then "
            "normalizes them to sum to one."
        ),
    )
    parser.add_argument(
        "--processes",
        "-p",
        type=int,
        default=None,
        help="Number of worker processes; default uses the system setting.",
    )
    parser.add_argument(
        "--residuals",
        "-r",
        action="store_true",
        help="Write per-sample L2 residual, RMSE, and number of fitted CpGs.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Write a stacked-bar plot of estimated fractions.",
    )
    parser.add_argument(
        "--slim",
        action="store_true",
        help="Write result values without row or column labels.",
    )
    parser.add_argument(
        "--out_prefix",
        "-o",
        type=Path,
        default=None,
        help=(
            "Output filename prefix. "
            "Default: use the input filename without extension."
        ),
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.processes is not None and args.processes < 1:
        parser.error("--processes must be at least 1")

    atlas_path = ATLAS_DIR / ATLAS_FILES[args.atlas]
    if not atlas_path.is_file():
        parser.error(
            f"Reference atlas was not found: {atlas_path}. "
            "Install or copy the reference files into CpGtools/src/deconv_data/."
        )

    try:
        Deconvolve(
            atlas_path=atlas_path,
            samples_path=args.samples_path,
            output_prefix=args.out_prefix,
            method=args.method,
            processes=args.processes,
            residuals=args.residuals,
            slim=args.slim,
            plot=args.plot,
        ).run()
    except DeconvolutionError as exc:
        parser.exit(1, f"Error: {exc}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

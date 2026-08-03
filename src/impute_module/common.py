# Portions of this file are derived from "fancyimpute", licensed under
# "Apache License, Version 2.0 "
#
# Modifications Copyright (c) 2024-2026 Liguo Wang
# Modified for use in CpGtools.
#
# See LICENSE.txt and the relevant third-party license notices for details.

from __future__ import annotations
import importlib
import logging
from collections.abc import Sequence
from typing import Any
import numpy as np
from numpy.typing import ArrayLike, NDArray

LOGGER = logging.getLogger(__name__)

FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]

def import_from(module_name: str, attribute_name: str) -> Any:
    """
    Import an attribute dynamically from a module.

    Parameters
    ----------
    module_name
        Fully qualified module name.

    attribute_name
        Name of the class, function, or object to import.

    Returns
    -------
    Any
        Imported attribute.

    Examples
    --------
    >>> GridSearchCV = import_from(
    ...     "sklearn.model_selection",
    ...     "GridSearchCV",
    ... )
    """
    if not isinstance(module_name, str) or not module_name.strip():
        raise ValueError("module_name must be a nonempty string.")

    if not isinstance(attribute_name, str) or not attribute_name.strip():
        raise ValueError("attribute_name must be a nonempty string.")

    try:
        imported_module = importlib.import_module(module_name)
    except ImportError as exc:
        raise ImportError(
            f"Could not import module {module_name!r}."
        ) from exc

    try:
        return getattr(imported_module, attribute_name)
    except AttributeError as exc:
        raise ImportError(
            f"Module {module_name!r} does not define "
            f"{attribute_name!r}."
        ) from exc


def _validate_masked_arrays(
    x_true: ArrayLike,
    x_pred: ArrayLike,
    mask: ArrayLike,
) -> tuple[FloatArray, FloatArray, BoolArray]:
    """Validate arrays used by masked error metrics."""
    true_values = np.asarray(x_true, dtype=float)
    predicted_values = np.asarray(x_pred, dtype=float)
    boolean_mask = np.asarray(mask)

    if true_values.shape != predicted_values.shape:
        raise ValueError(
            "x_true and x_pred must have the same shape; "
            f"received {true_values.shape} and {predicted_values.shape}."
        )

    if boolean_mask.shape != true_values.shape:
        raise ValueError(
            "mask must have the same shape as x_true and x_pred; "
            f"received {boolean_mask.shape} and {true_values.shape}."
        )

    if boolean_mask.dtype != np.bool_:
        raise TypeError("mask must contain boolean values.")

    if not boolean_mask.any():
        raise ValueError("mask does not select any values.")

    true_selected = true_values[boolean_mask]
    predicted_selected = predicted_values[boolean_mask]

    if not np.isfinite(true_selected).all():
        raise ValueError(
            "x_true contains non-finite values at selected positions."
        )

    if not np.isfinite(predicted_selected).all():
        raise ValueError(
            "x_pred contains non-finite values at selected positions."
        )

    return true_values, predicted_values, boolean_mask


def masked_mae(
    x_true: ArrayLike,
    x_pred: ArrayLike,
    mask: ArrayLike,
) -> float:
    """
    Calculate mean absolute error over positions selected by ``mask``.
    """
    true_values, predicted_values, boolean_mask = _validate_masked_arrays(
        x_true,
        x_pred,
        mask,
    )

    return float(
        np.mean(
            np.abs(
                true_values[boolean_mask]
                - predicted_values[boolean_mask]
            )
        )
    )


def masked_mse(
    x_true: ArrayLike,
    x_pred: ArrayLike,
    mask: ArrayLike,
) -> float:
    """
    Calculate mean squared error over positions selected by ``mask``.
    """
    true_values, predicted_values, boolean_mask = _validate_masked_arrays(
        x_true,
        x_pred,
        mask,
    )

    differences = (
        true_values[boolean_mask]
        - predicted_values[boolean_mask]
    )

    return float(np.mean(np.square(differences)))


def generate_random_column_samples(
    column: ArrayLike,
    *,
    random_state: int | np.random.Generator | None = None,
) -> FloatArray:
    """
    Generate replacement values for missing entries in one column.

    Samples are drawn from a normal distribution estimated from the observed
    values. If the observed values have zero variance, their mean is returned
    for every missing position.

    Parameters
    ----------
    column
        One-dimensional numeric array containing optional missing values.

    random_state
        Integer seed, NumPy generator, or ``None``.

    Returns
    -------
    numpy.ndarray
        One-dimensional array containing one generated value per missing
        entry.

    Raises
    ------
    ValueError
        If the input is not one-dimensional or contains no observed values.
    """
    values = np.asarray(column, dtype=float)

    if values.ndim != 1:
        raise ValueError("column must be one-dimensional.")

    if np.isinf(values).any():
        raise ValueError("column must not contain infinite values.")

    missing_mask = np.isnan(values)
    number_missing = int(missing_mask.sum())

    if number_missing == 0:
        return np.empty(0, dtype=float)

    observed = values[~missing_mask]

    if observed.size == 0:
        raise ValueError(
            "Cannot generate samples because the column contains "
            "no observed values."
        )

    mean = float(observed.mean())
    standard_deviation = float(observed.std(ddof=0))

    if np.isclose(standard_deviation, 0.0):
        return np.full(number_missing, mean, dtype=float)

    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    return rng.normal(
        loc=mean,
        scale=standard_deviation,
        size=number_missing,
    )


def choose_solution_using_percentiles(
    x_original: ArrayLike,
    solutions: Sequence[ArrayLike],
    *,
    parameters: Sequence[Any] | None = None,
    percentiles: Sequence[float] = tuple(range(10, 100, 10)),
    verbose: bool = False,
) -> FloatArray:
    """
    Select the candidate whose imputed-value percentiles most closely match
    the observed-value percentiles.

    For every candidate and every eligible column, this function compares
    the requested percentiles of the originally missing entries with those
    of the originally observed entries. The candidate score is the mean of
    these column-level mean squared errors.

    Parameters
    ----------
    x_original
        Original two-dimensional matrix containing missing values.

    solutions
        Candidate completed matrices. Every candidate must have the same
        shape as ``x_original`` and contain only finite values.

    parameters
        Optional parameter labels corresponding to the candidate solutions.
        These are used only in progress messages.

    percentiles
        Percentiles to compare. Values must lie between 0 and 100.

    verbose
        Emit one informational log message per candidate.

    Returns
    -------
    numpy.ndarray
        Copy of the candidate matrix with the lowest average percentile MSE.

    Raises
    ------
    ValueError
        If the inputs are invalid or no candidate can be evaluated.
    """
    original = np.asarray(x_original, dtype=float)

    if original.ndim != 2:
        raise ValueError("x_original must be two-dimensional.")

    if original.size == 0:
        raise ValueError("x_original must not be empty.")

    if np.isinf(original).any():
        raise ValueError("x_original must not contain infinite values.")

    if not solutions:
        raise ValueError("solutions must contain at least one candidate.")

    if parameters is not None and len(parameters) != len(solutions):
        raise ValueError(
            "parameters must have the same length as solutions."
        )

    percentile_values = np.asarray(percentiles, dtype=float)

    if percentile_values.ndim != 1 or percentile_values.size == 0:
        raise ValueError(
            "percentiles must be a nonempty one-dimensional sequence."
        )

    if not np.isfinite(percentile_values).all():
        raise ValueError("percentiles must contain only finite values.")

    if ((percentile_values < 0) | (percentile_values > 100)).any():
        raise ValueError(
            "percentiles must contain values between 0 and 100."
        )

    missing_mask = np.isnan(original)

    if not missing_mask.any():
        raise ValueError(
            "x_original does not contain missing values, so candidate "
            "solutions cannot be compared at missing positions."
        )

    best_solution: FloatArray | None = None
    best_score = np.inf

    for candidate_index, candidate in enumerate(solutions):
        candidate_values = np.asarray(candidate, dtype=float)

        if candidate_values.shape != original.shape:
            raise ValueError(
                f"Candidate {candidate_index + 1} has shape "
                f"{candidate_values.shape}; expected {original.shape}."
            )

        if not np.isfinite(candidate_values).all():
            raise ValueError(
                f"Candidate {candidate_index + 1} contains missing or "
                "infinite values."
            )

        column_errors: list[float] = []

        for column_index in range(original.shape[1]):
            column_missing = missing_mask[:, column_index]
            column_observed = ~column_missing

            # Percentile comparisons are unstable with fewer than two values.
            if column_missing.sum() < 2 or column_observed.sum() < 2:
                continue

            missing_values = candidate_values[
                column_missing,
                column_index,
            ]
            observed_values = candidate_values[
                column_observed,
                column_index,
            ]

            missing_percentiles = np.percentile(
                missing_values,
                percentile_values,
            )
            observed_percentiles = np.percentile(
                observed_values,
                percentile_values,
            )

            column_mse = float(
                np.mean(
                    np.square(
                        missing_percentiles - observed_percentiles
                    )
                )
            )
            column_errors.append(column_mse)

        if not column_errors:
            candidate_score = np.inf
        else:
            candidate_score = float(np.mean(column_errors))

        if verbose:
            parameter_text = (
                f", parameter={parameters[candidate_index]!r}"
                if parameters is not None
                else ""
            )

            LOGGER.info(
                "Candidate %d/%d%s: percentile MSE=%.6g, "
                "eligible columns=%d",
                candidate_index + 1,
                len(solutions),
                parameter_text,
                candidate_score,
                len(column_errors),
            )

        if candidate_score < best_score:
            best_score = candidate_score
            best_solution = candidate_values.copy()

    if best_solution is None:
        raise ValueError(
            "No candidate could be evaluated. Each eligible column needs "
            "at least two originally missing values and two observed values."
        )

    if verbose:
        LOGGER.info(
            "Selected candidate with percentile MSE %.6g.",
            best_score,
        )

    return best_solution

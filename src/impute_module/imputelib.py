# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: wangliguo78@gmail.com
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

from __future__ import annotations
import logging
import pandas as pd
import numpy as np
from typing import Literal
from pathlib import Path
from scipy.spatial import KDTree
from impute_module.softimpute import SoftImpute
from sklearn.impute import KNNImputer
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from numbers import Integral

from collections.abc import Callable

LOGGER = logging.getLogger("nafiller")
Statistic = Literal["mean", "median", "min", "max"]
Axis = Literal["index", "columns"]
Weights = Literal["uniform", "distance"]
SearchAxis = Literal["index", "columns"]
STATISTICS = {
    "mean": pd.DataFrame.mean,
    "median": pd.DataFrame.median,
    "min": pd.DataFrame.min,
    "max": pd.DataFrame.max,
}
SUMMARY_FUNCTIONS = {
    "mean": np.mean,
    "median": np.median,
}

def configure_logging(debug: bool = False) -> None:
    """Configure command-line logging."""
    logging.basicConfig(
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        level=logging.DEBUG if debug else logging.INFO,
    )

def _nan_indices(data):
    """ Finds the indices of all missing values.

    Parameters
    ----------
    data: numpy.ndarray

    Returns
    -------
    List of tuples
        Indices of all missing values in tuple format; (i, j)
    """
    return np.argwhere(np.isnan(data))

def _moving_window(
    data: np.ndarray | pd.DataFrame,
    nindex: int | None = None,
    wsize: int = 5,
    errors: str = "coerce",
    func: Callable = np.mean,
    inplace: bool = False,) -> np.ndarray | pd.DataFrame:
    """
    Impute missing values using neighboring values in a horizontal window.

    The function processes each row independently. For every missing value,
    it extracts a window from the same row and replaces the missing value
    with a statistic calculated from the observed values in that window.

    To apply the window vertically to DataFrame columns, transpose the input
    before calling this function and transpose the result afterward.

    Parameters
    ----------
    data
        Two-dimensional NumPy array or pandas DataFrame.

    nindex
        Position of the missing value within the window.

        - ``None``: center the missing value in the window.
        - ``0``: use values to the right.
        - ``-1``: use values to the left.

    wsize
        Total size of the moving window, including the missing-value
        position. When ``nindex=None``, ``wsize`` must be odd.

    errors
        How to handle a missing value for which no usable neighboring values
        are found:

        - ``"raise"``: raise an exception.
        - ``"coerce"``: retry with a centered window.
        - ``"ignore"``: leave the value missing.

    func
        Callable used to summarize observed values in the window, such as
        ``numpy.mean`` or ``numpy.median``.

    inplace
        If ``True`` and the input is a NumPy array, modify it directly.
        Otherwise, operate on a copy.

    Returns
    -------
    numpy.ndarray or pandas.DataFrame
        Imputed data in the same type and shape as the input.

    Raises
    ------
    TypeError
        If the input or summary function has an invalid type.

    ValueError
        If an argument is invalid or imputation fails with
        ``errors="raise"``.
    """
    if not isinstance(data, (np.ndarray, pd.DataFrame)):
        raise TypeError("data must be a NumPy array or pandas DataFrame.")

    if not callable(func):
        raise TypeError("func must be callable.")

    if errors not in ("raise", "coerce", "ignore"):
        raise ValueError(
            "errors must be one of 'raise', 'coerce', or 'ignore'."
        )

    if not isinstance(wsize, int) or isinstance(wsize, bool) or wsize < 1:
        raise ValueError("wsize must be a positive integer.")

    if nindex is not None:
        if not isinstance(nindex, int) or isinstance(nindex, bool):
            raise TypeError("nindex must be an integer or None.")

        if nindex != -1 and not 0 <= nindex < wsize:
            raise ValueError(
                "nindex must be -1 or an integer from 0 to wsize - 1."
            )

    if nindex is None and wsize % 2 == 0:
        raise ValueError(
            "wsize must be odd when nindex is None because the missing "
            "value is placed at the center."
        )

    is_dataframe = isinstance(data, pd.DataFrame)

    if is_dataframe:
        index = data.index
        columns = data.columns
        values = data.to_numpy(dtype=float, copy=True)
    else:
        values = np.asarray(data, dtype=float)

        if values.ndim != 2:
            raise ValueError("data must be two-dimensional.")

        if not inplace:
            values = values.copy()

    if values.ndim != 2:
        raise ValueError("data must be two-dimensional.")

    if nindex is None:
        left_width = wsize // 2
        right_width = wsize // 2
    elif nindex == -1:
        left_width = wsize - 1
        right_width = 0
    else:
        left_width = nindex
        right_width = wsize - nindex - 1

    n_columns = values.shape[1]

    while True:
        missing_locations = np.argwhere(np.isnan(values))

        if missing_locations.size == 0:
            break

        missing_before = len(missing_locations)

        for row_index, column_index in missing_locations:
            left = max(0, column_index - left_width)
            right = min(
                n_columns,
                column_index + right_width + 1,
            )

            window = values[row_index, left:right]
            observed = window[~np.isnan(window)]

            if observed.size > 0:
                try:
                    estimate = func(observed)

                    if np.isfinite(estimate):
                        values[row_index, column_index] = estimate
                        continue
                except Exception:
                    if errors == "raise":
                        raise

            if errors == "coerce":
                fallback_width = wsize // 2

                fallback_left = max(
                    0,
                    column_index - fallback_width,
                )
                fallback_right = min(
                    n_columns,
                    column_index + fallback_width + 1,
                )

                fallback_window = values[
                    row_index,
                    fallback_left:fallback_right,
                ]
                fallback_observed = fallback_window[
                    ~np.isnan(fallback_window)
                ]

                if fallback_observed.size > 0:
                    try:
                        estimate = func(fallback_observed)

                        if np.isfinite(estimate):
                            values[row_index, column_index] = estimate
                            continue
                    except Exception:
                        if errors == "raise":
                            raise

            if errors == "raise":
                raise ValueError(
                    "Unable to impute the missing value at "
                    f"row {row_index}, column {column_index}: "
                    "the moving window contains no usable values."
                )

            # errors == "ignore": leave the value unchanged.

        missing_after = int(np.isnan(values).sum())

        # Stop when no missing value was filled during this iteration.
        if missing_after == missing_before:
            break

    if is_dataframe:
        return pd.DataFrame(
            values,
            index=index,
            columns=columns,
        )

    return values


def _buck_iterative(
    data: np.ndarray,
    eps: float = 1e-3,
    max_iter: int = 200,
    random_state: int | None = None,
    copy: bool = True,) -> np.ndarray:
    """
    Impute missing values using iterative Buck regression imputation.

    Missing values are first initialized with their row means. Each
    column that originally contained missing values is then regressed on
    all other columns. The fitted model predicts the originally missing
    entries in that column. This process is repeated until the imputed
    values converge or ``max_iter`` is reached.

    Parameters
    ----------
    data
        Two-dimensional numeric array containing missing values represented
        by ``numpy.nan``. Rows are observations and columns are variables.

    eps
        Convergence tolerance. Iteration stops when the largest relative
        change among all imputed values is less than ``eps``.

        Relative change is calculated as::

            abs(current - previous) / max(abs(previous), 1e-12)

    max_iter
        Maximum number of full iterations over columns containing missing
        values.

    random_state
        Random seed used to shuffle the order in which incomplete columns
        are updated during each iteration. If ``None``, columns are processed
        in their original order.

    copy
        If ``True``, operate on a copy and leave the input unchanged.
        If ``False``, the supplied floating-point array may be modified.

    Returns
    -------
    numpy.ndarray
        Array with missing values imputed.

    Raises
    ------
    TypeError
        If an argument has an invalid type.

    ValueError
        If the input is invalid, a column is entirely missing, or regression
        cannot be performed for an incomplete column.

    RuntimeError
        If the algorithm does not converge within ``max_iter`` iterations.

    Notes
    -----
    This is an iterative regression implementation based on:

    Buck, S. F. (1960). A Method of Estimation of Missing Values in
    Multivariate Data Suitable for Use with an Electronic Computer.
    Journal of the Royal Statistical Society: Series B, 22(2), 302-306.
    """
    values = np.asarray(data)

    if values.ndim != 2:
        raise ValueError("data must be a two-dimensional array.")

    if not np.issubdtype(values.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if not isinstance(eps, (int, float)) or isinstance(eps, bool):
        raise TypeError("eps must be numeric.")

    if eps <= 0:
        raise ValueError("eps must be greater than zero.")

    if (
        not isinstance(max_iter, (int, np.integer))
        or isinstance(max_iter, bool)):
        raise TypeError("max_iter must be an integer.")

    if max_iter < 1:
        raise ValueError("max_iter must be a positive integer.")

    # Ensure NaN values can be represented.
    if copy or not np.issubdtype(values.dtype, np.floating):
        values = values.astype(float, copy=True)

    if np.isinf(values).any():
        raise ValueError(
            "data contains infinite values; replace or remove them first.")

    missing_mask = np.isnan(values)

    if not missing_mask.any():
        return values

    all_missing_columns = missing_mask.all(axis=0)

    if all_missing_columns.any():
        columns = np.flatnonzero(all_missing_columns).tolist()
        raise ValueError(
            "Buck imputation cannot estimate columns containing only "
            f"missing values: {columns}.")

    incomplete_columns = np.flatnonzero(missing_mask.any(axis=0))

    if values.shape[1] < 2:
        raise ValueError(
            "Buck regression imputation requires at least two columns.")

    # Initial imputation with column means.
    #column_means = np.nanmean(values, axis=0)
    #missing_rows, missing_columns = np.where(missing_mask)
    #values[missing_rows, missing_columns] = column_means[missing_columns]


    # Initial imputation with row means.
    row_means = np.nanmean(values, axis=1)
    all_missing_rows = np.isnan(row_means)
    if all_missing_rows.any():
        rows = np.flatnonzero(all_missing_rows).tolist()
        raise ValueError(
            "Row-mean initialization cannot be performed for rows containing "
            f"only missing values: {rows}."
        )
    missing_rows, missing_columns = np.where(missing_mask)
    values[missing_rows, missing_columns] = row_means[missing_rows]


    rng = (
        np.random.default_rng(random_state)
        if random_state is not None
        else None)

    for iteration in range(1, max_iter + 1):
        previous = values[missing_mask].copy()

        column_order = incomplete_columns.copy()
        if rng is not None:
            rng.shuffle(column_order)

        for dependent_column in column_order:
            originally_missing = missing_mask[:, dependent_column]
            training_rows = ~originally_missing

            if training_rows.sum() < 2:
                raise ValueError(
                    f"Column {dependent_column} has fewer than two observed "
                    "values and cannot be fitted reliably.")

            predictor_columns = np.arange(values.shape[1]) != dependent_column

            x_train = values[np.ix_(training_rows, predictor_columns)]
            y_train = values[training_rows,dependent_column]

            x_predict = values[np.ix_(originally_missing, predictor_columns)]

            if x_predict.shape[0] == 0:
                continue

            model = LinearRegression()
            model.fit(x_train, y_train)

            values[originally_missing,dependent_column] = model.predict(x_predict)

        current = values[missing_mask]

        denominator = np.maximum(np.abs(previous), 1e-12)
        relative_change = np.abs(current - previous) / denominator
        max_change = float(np.max(relative_change))

        if max_change < eps:
            return values

    raise RuntimeError(
        "Buck imputation did not converge after "
        f"{max_iter} iterations. Consider increasing max_iter or eps.")
    
def _rf_iterative(
    data: np.ndarray,
    *,
    n_estimators: int = 100,
    max_iter: int = 10,
    eps: float = 1e-3,
    max_depth: int | None = None,
    min_samples_leaf: int = 1,
    random_state: int | None = 1234,
    n_jobs: int | None = -1,
    min_value: float = -np.inf,
    max_value: float = np.inf,) -> np.ndarray:
    """
    Impute missing values using iterative random-forest regression.

    Missing values are initially replaced with the mean of their row. Each
    column containing missing values is then modeled as a regression target
    using the remaining columns as predictors. Originally missing values are
    repeatedly updated until convergence or until ``max_iter`` is reached.

    Parameters
    ----------
    data
        Two-dimensional numeric array. Rows are observations and columns are
        regression variables.

    n_estimators
        Number of trees fitted for each incomplete column.

    max_iter
        Maximum number of complete imputation rounds.

    eps
        Convergence threshold. Iteration stops when the maximum normalized
        change among imputed values is less than ``eps``.

    max_features
        Number or fraction of predictor variables considered at each tree
        split.

    max_depth
        Maximum depth of each tree. ``None`` permits unrestricted depth.

    min_samples_leaf
        Minimum number of training observations in each terminal leaf.

    random_state
        Base random seed. Use an integer for reproducible results.

    n_jobs
        Number of parallel jobs used to fit each random forest. ``-1`` uses
        all available processors.

    min_value
        Lower bound applied to imputed values.

    max_value
        Upper bound applied to imputed values.

    Returns
    -------
    numpy.ndarray
        A new array containing the imputed values.

    Raises
    ------
    TypeError
        If ``data`` is not a NumPy array.

    ValueError
        If the input or parameters are invalid, or if imputation cannot be
        performed.
    """
    if not isinstance(data, np.ndarray):
        raise TypeError("data must be a NumPy array.")

    if data.ndim != 2:
        raise ValueError("data must be a two-dimensional array.")

    if data.size == 0:
        raise ValueError("data must not be empty.")

    if not np.issubdtype(data.dtype, np.number):
        raise TypeError("data must contain numeric values.")

    if np.isinf(data).any():
        raise ValueError("data must not contain infinite values.")

    if (
        not isinstance(n_estimators, (int, np.integer))
        or isinstance(n_estimators, bool)
        or n_estimators < 1
    ):
        raise ValueError("n_estimators must be a positive integer.")

    if (
        not isinstance(max_iter, (int, np.integer))
        or isinstance(max_iter, bool)
        or max_iter < 1
    ):
        raise ValueError("max_iter must be a positive integer.")

    if (
        not isinstance(eps, (int, float, np.integer, np.floating))
        or isinstance(eps, bool)
        or not np.isfinite(eps)
        or eps <= 0
    ):
        raise ValueError("eps must be a finite number greater than zero.")

    if (
        not isinstance(min_samples_leaf, (int, np.integer))
        or isinstance(min_samples_leaf, bool)
        or min_samples_leaf < 1
    ):
        raise ValueError("min_samples_leaf must be a positive integer.")

    if min_value >= max_value:
        raise ValueError("min_value must be smaller than max_value.")

    values = np.array(data, dtype=float, copy=True)
    missing_mask = np.isnan(values)

    if not missing_mask.any():
        return values

    n_rows, n_columns = values.shape

    if n_columns < 2:
        raise ValueError(
            "Random-forest imputation requires at least two columns.")

    # Every row needs at least one observed value for row-mean initialization.
    all_missing_rows = missing_mask.all(axis=1)
    if all_missing_rows.any():
        row_numbers = np.flatnonzero(all_missing_rows).tolist()
        raise ValueError(
            "Cannot initialize rows containing only missing values. "
            f"Zero-based row positions: {row_numbers}.")

    # A completely missing target column has no observed response values from
    # which a random-forest model can be trained.
    all_missing_columns = missing_mask.all(axis=0)
    if all_missing_columns.any():
        column_numbers = np.flatnonzero(all_missing_columns).tolist()
        raise ValueError(
            "Cannot model columns containing only missing values. "
            f"Zero-based column positions: {column_numbers}.")

    # Initial imputation with row means.
    row_means = np.nanmean(values, axis=1)
    missing_rows, missing_columns = np.where(missing_mask)
    values[missing_rows, missing_columns] = row_means[missing_rows]

    # Only columns that originally contained missing values need models.
    incomplete_columns = np.flatnonzero(missing_mask.any(axis=0))

    # Start with columns containing fewer missing entries. Their models have
    # more observed response values available for training.
    missing_counts = missing_mask[:, incomplete_columns].sum(axis=0)
    incomplete_columns = incomplete_columns[np.argsort(missing_counts)]

    for iteration in range(1, max_iter + 1):
        previous_missing_values = values[missing_mask].copy()

        for target_column in incomplete_columns:
            target_missing = missing_mask[:, target_column]
            target_observed = ~target_missing

            predictor_columns = np.arange(n_columns) != target_column

            x_train = values[target_observed][:, predictor_columns]
            y_train = values[target_observed, target_column]
            x_predict = values[target_missing][:, predictor_columns]

            if y_train.size < 2:
                raise ValueError(
                    "Column "
                    f"{target_column} has fewer than two observed values and "
                    "cannot be modeled reliably.")

            if x_train.shape[1] == 0:
                raise ValueError(
                    f"No predictors are available for column {target_column}.")

            # Give every fitted forest a deterministic but distinct seed.
            model_seed = (
                None
                if random_state is None
                else int(random_state + iteration * n_columns + target_column))

            model = RandomForestRegressor(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_leaf=min_samples_leaf,
                random_state=model_seed,
                n_jobs=n_jobs,)

            model.fit(x_train, y_train)
            predictions = model.predict(x_predict)

            lower = -np.inf if min_value is None else min_value
            upper = np.inf if max_value is None else max_value
            predictions = np.clip(predictions, lower, upper)

            values[target_missing, target_column] = predictions

        current_missing_values = values[missing_mask]

        absolute_change = np.abs(
            current_missing_values - previous_missing_values)

        scale = np.maximum(
            np.abs(previous_missing_values),
            np.finfo(float).eps,)

        max_relative_change = float(
            np.max(absolute_change / scale))

        LOGGER.debug(
            "Random-forest iteration %d: maximum relative change = %.6g",
            iteration,
            max_relative_change,)

        if max_relative_change < eps:
            LOGGER.info(
                "Random-forest imputation converged after %d iteration(s).",
                iteration,)
            break
    else:
        LOGGER.warning(
            "Random-forest imputation did not converge after %d iterations.",
            max_iter,)

    # Guard against accidentally modifying observed entries.
    values[~missing_mask] = data[~missing_mask]

    if np.isnan(values).any():
        raise RuntimeError(
            "Random-forest imputation finished with missing values remaining.")

    return values


def _external_ref(
    data: np.ndarray,
    na_locations: np.ndarray,
    ref_data: np.ndarray,
    k: int = 3,
    eps: float = 0.0,
    p: float = 2.0,
    distance_upper_bound: float = np.inf,
    leafsize: int = 10,) -> np.ndarray:
    """
    Impute missing values using nearest neighbors from an external reference.

    Nearest neighbors are identified using only features that were originally
    observed in each query row. Features listed in ``na_locations`` are
    excluded from the distance calculation.

    Missing values are replaced with an inverse-distance-weighted average of
    the corresponding values from the selected reference neighbors.

    Parameters
    ----------
    data
        Finite query matrix with observations in rows and features in
        columns. Originally missing values may have been temporarily filled
        before calling this function.

    na_locations
        Locations of the originally missing values. Accepted forms are:

        - an ``(n_missing, 2)`` array of ``(row, column)`` coordinates;
        - the tuple returned by ``numpy.where``;
        - the tuple returned by ``numpy.nonzero``.

    ref_data
        Complete external reference matrix. Reference observations must be
        in rows and features in columns. It must not contain NaN or infinite
        values.

    k
        Number of nearest reference neighbors.

    eps
        Approximation factor passed to ``scipy.spatial.KDTree.query``.
        Zero performs an exact nearest-neighbor search.

    p
        Minkowski distance norm:

        - ``1``: Manhattan distance
        - ``2``: Euclidean distance
        - ``numpy.inf``: Chebyshev distance

    distance_upper_bound
        Maximum allowed neighbor distance.

    leafsize
        KDTree leaf size.

    Returns
    -------
    numpy.ndarray
        A copy of ``data`` with values at ``na_locations`` replaced by
        reference-based KNN estimates.

    Raises
    ------
    TypeError
        If an integer argument has an invalid type.

    ValueError
        If matrix dimensions or parameter values are invalid, reference data
        are nonfinite, or no usable neighbors can be found.
    """
    data = np.asarray(data, dtype=float)
    ref_data = np.asarray(ref_data, dtype=float)

    if data.ndim != 2:
        raise ValueError("data must be a two-dimensional array.")

    if ref_data.ndim != 2:
        raise ValueError("ref_data must be a two-dimensional array.")

    if data.shape[1] != ref_data.shape[1]:
        raise ValueError(
            "data and ref_data must contain the same number of features."
        )

    if not np.isfinite(data).all():
        raise ValueError(
            "data must be finite. Temporarily fill missing query values "
            "before calling external_ref()."
        )

    if not np.isfinite(ref_data).all():
        raise ValueError(
            "ref_data must not contain NaN or infinite values."
        )

    if not isinstance(k, (int, np.integer)) or isinstance(k, bool):
        raise TypeError("k must be an integer.")

    if k < 1:
        raise ValueError("k must be greater than or equal to 1.")

    if k > ref_data.shape[0]:
        raise ValueError(
            f"k={k} exceeds the number of reference observations "
            f"({ref_data.shape[0]}).")

    if eps < 0:
        raise ValueError("eps must be nonnegative.")

    if p < 1:
        raise ValueError("p must be greater than or equal to 1.")

    if distance_upper_bound <= 0:
        raise ValueError(
            "distance_upper_bound must be positive or numpy.inf.")

    if (
        not isinstance(leafsize, (int, np.integer))
        or isinstance(leafsize, bool)):
        raise TypeError("leafsize must be an integer.")

    if leafsize < 1:
        raise ValueError("leafsize must be greater than or equal to 1.")

    na_coordinates = _normalize_na_locations(
        na_locations=na_locations,
        shape=data.shape,)

    output = data.copy()

    if na_coordinates.size == 0:
        return output

    # Process every query observation only once. All originally missing
    # features in the same observation use the same neighbor set.
    for row_index in np.unique(na_coordinates[:, 0]):
        missing_columns = na_coordinates[
            na_coordinates[:, 0] == row_index, 1]

        missing_mask = np.zeros(data.shape[1], dtype=bool)
        missing_mask[missing_columns] = True
        observed_mask = ~missing_mask

        if not observed_mask.any():
            raise ValueError(
                f"Query observation {row_index} has no observed features "
                "available for calculating distances.")

        query_features = data[row_index, observed_mask]
        reference_features = ref_data[:, observed_mask]

        # A separate tree is necessary because each query observation may
        # have a different set of originally missing features.
        tree = KDTree(
            reference_features,
            leafsize=leafsize,)

        distances, indices = tree.query(
            query_features,
            k=k,
            eps=eps,
            p=p,
            distance_upper_bound=distance_upper_bound,)

        # KDTree returns scalars when k == 1.
        distances = np.atleast_1d(distances).astype(float)
        indices = np.atleast_1d(indices).astype(int)

        # When distance_upper_bound excludes a neighbor, KDTree represents
        # it with an infinite distance and an out-of-range index.
        valid = (
            np.isfinite(distances)
            & (indices >= 0)
            & (indices < ref_data.shape[0]))

        distances = distances[valid]
        indices = indices[valid]

        if indices.size == 0:
            raise ValueError(
                f"No reference neighbors were found for query observation "
                f"{row_index} within distance_upper_bound="
                f"{distance_upper_bound}.")

        weights = _inverse_distance_weights(distances)

        neighbor_values = ref_data[
            np.ix_(indices, missing_columns)]

        output[row_index, missing_columns] = (
            weights @ neighbor_values)

    return output

def _random_impute(
    data: np.ndarray,
    *,
    random_state: int | np.random.Generator | None = None,
    inplace: bool = False,) -> np.ndarray:
    """
    Impute missing values by randomly sampling observed values from the
    same column.

    Each missing value is replaced by a randomly selected observed value
    from its own column. Sampling is performed with replacement.

    Parameters
    ----------
    data
        Two-dimensional numeric array.

    random_state
        Integer seed, NumPy Generator, or ``None``.

    inplace
        If ``True``, modify the supplied floating-point array directly.
        Otherwise operate on a copy.

    Returns
    -------
    numpy.ndarray
        Imputed array.

    Raises
    ------
    TypeError
        If ``data`` is not a NumPy array.

    ValueError
        If the input is invalid or a column contains only missing values.
    """
    if not isinstance(data, np.ndarray):
        raise TypeError("data must be a NumPy array.")
    if data.ndim != 2:
        raise ValueError("data must be two-dimensional.")
    if data.size == 0:
        raise ValueError("data must not be empty.")
    if np.isinf(data).any():
        raise ValueError("data must not contain infinite values.")
    if inplace:
        if not np.issubdtype(data.dtype, np.floating):
            raise TypeError(
                "inplace=True requires a floating-point array.")
        values = data
    else:
        values = np.array(data, dtype=float, copy=True)

    missing_mask = np.isnan(values)

    if not missing_mask.any():
        return values

    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    for column in range(values.shape[1]):
        missing = missing_mask[:, column]
        if not missing.any():
            continue
        observed = values[~missing, column]
        if observed.size == 0:
            raise ValueError(
                f"Column {column} contains only missing values.")
        values[missing, column] = rng.choice(
            observed,
            size=int(missing.sum()),
            replace=True,)

    return values

def _normalize_na_locations(
    na_locations,
    shape: tuple[int, int],) -> np.ndarray:
    """Convert supported missing-location formats to an (n, 2) array."""
    if isinstance(na_locations, tuple):
        if len(na_locations) != 2:
            raise ValueError(
                "A tuple of na_locations must contain row and column arrays."
            )

        rows = np.asarray(na_locations[0], dtype=int)
        columns = np.asarray(na_locations[1], dtype=int)

        if rows.shape != columns.shape:
            raise ValueError(
                "The row and column arrays in na_locations must have "
                "matching shapes.")

        coordinates = np.column_stack((rows, columns))
    else:
        coordinates = np.asarray(na_locations, dtype=int)

        if coordinates.size == 0:
            return np.empty((0, 2), dtype=int)

        if coordinates.ndim != 2 or coordinates.shape[1] != 2:
            raise ValueError(
                "na_locations must be an (n_missing, 2) array or a "
                "(rows, columns) tuple.")

    if coordinates.size == 0:
        return np.empty((0, 2), dtype=int)

    if (
        (coordinates[:, 0] < 0).any()
        or (coordinates[:, 0] >= shape[0]).any()
        or (coordinates[:, 1] < 0).any()
        or (coordinates[:, 1] >= shape[1]).any()):
        raise ValueError(
            "na_locations contains coordinates outside the data matrix.")

    return np.unique(coordinates, axis=0)


def _inverse_distance_weights(distances: np.ndarray,) -> np.ndarray:
    """
    Calculate normalized inverse-distance weights.

    If one or more reference observations are exact matches, only those
    exact matches receive weight.
    """
    distances = np.asarray(distances, dtype=float)

    exact_matches = np.isclose(
        distances,
        0.0,
        rtol=0.0,
        atol=np.finfo(float).eps,
    )

    if exact_matches.any():
        weights = exact_matches.astype(float)
    else:
        weights = 1.0 / distances

    return weights / weights.sum()


def _prepare_knn_query(
    input_array: np.ndarray,
    ref_array: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Create a finite query matrix for KNN while preserving NA locations.

    Missing query values are temporarily filled with feature means from the
    input matrix. If a feature is entirely missing, its reference mean is
    used instead.
    """
    input_array = np.asarray(input_array, dtype=float)
    ref_array = np.asarray(ref_array, dtype=float)

    if np.isinf(input_array).any():
        raise ValueError(
            "The input matrix contains infinite values. "
            "Replace or remove them before KNN imputation."
        )

    if not np.isfinite(ref_array).all():
        raise ValueError(
            "The external reference matrix must contain only finite values."
        )

    na_locations = _nan_indices(input_array)
    query_array = input_array.copy()

    # Calculate one temporary replacement value per feature.
    with np.errstate(all="ignore"):
        feature_means = np.nanmean(query_array, axis=0)

    # Use reference means for features that are entirely missing in input.
    missing_feature_means = np.isnan(feature_means)
    if missing_feature_means.any():
        reference_means = np.mean(ref_array, axis=0)
        feature_means[missing_feature_means] = reference_means[
            missing_feature_means
        ]

    if not np.isfinite(feature_means).all():
        raise ValueError(
            "A finite temporary value could not be calculated for one or "
            "more features."
        )

    rows, columns = np.where(np.isnan(query_array))
    query_array[rows, columns] = feature_means[columns]

    return query_array, na_locations

def report_missing_values(
    df: pd.DataFrame,
    row_output: str | Path,
    column_output: str | Path,) -> None:
    """
    Count missing values and write row- and column-level reports.

    Parameters
    ----------
    df
        Input matrix with CpGs in rows and samples in columns.

    row_output
        Output file for missing-value counts per row. The file contains
        two tab-delimited columns: ``CpG`` and ``missing_count``.

    column_output
        Output file for missing-value counts per column. The file contains
        two tab-delimited columns: ``sample`` and ``missing_count``.

    Returns
    -------
    int
        Total number of missing values in the matrix.

    Raises
    ------
    TypeError
        If ``df`` is not a pandas DataFrame.

    ValueError
        If the row or column labels contain duplicates.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if df.index.has_duplicates:
        raise ValueError("The input matrix contains duplicated row labels.")

    if df.columns.has_duplicates:
        raise ValueError("The input matrix contains duplicated column labels.")

    row_output = Path(row_output)
    column_output = Path(column_output)

    row_output.parent.mkdir(parents=True, exist_ok=True)
    column_output.parent.mkdir(parents=True, exist_ok=True)

    missing_mask = df.isna()

    missing_by_row = missing_mask.sum(axis="columns").astype(int)
    missing_by_column = missing_mask.sum(axis="index").astype(int)

    row_report = missing_by_row.rename("missing_count").rename_axis("CpG")
    column_report = (
        missing_by_column.rename("missing_count").rename_axis("sample")
    )

    row_report.to_csv(row_output, sep="\t", header=True)
    column_report.to_csv(column_output, sep="\t", header=True)
    

def fill_ref(
    df: pd.DataFrame,
    ref: pd.DataFrame,
    axis: SearchAxis = "columns",
    k: int = 3,
    eps: float = 0.0,
    p: float = 2.0,
    distance_upper_bound: float = np.inf,
    leafsize: int = 10,) -> pd.DataFrame:
    """
    Impute missing values using nearest neighbors from an external reference.

    This method is similar to fast KNN imputation, except that nearest
    neighbors are searched in an external reference matrix rather than in
    the input matrix. Missing values are replaced using the corresponding
    values from the selected reference neighbors.

    Parameters
    ----------
    df
        Input matrix containing missing values. For DNA methylation data,
        CpGs are expected in rows and samples in columns.

    ref
        Complete external reference matrix. It must not contain missing
        values.

        When ``axis="columns"``, columns are treated as observations and
        rows as features. The input and reference matrices are aligned by
        row labels.

        When ``axis="index"``, rows are treated as observations and columns
        as features. The input and reference matrices are aligned by column
        labels.

    axis
        Dimension in which nearest neighbors are searched:

        - ``"columns"``: search reference columns, typically samples.
          This is the default.
        - ``"index"``: search reference rows, typically CpGs.

    k
        Number of nearest neighbors used for imputation. Must be a positive
        integer and cannot exceed the number of candidate observations in
        the reference matrix.

    eps
        Approximation factor passed to ``scipy.spatial.KDTree.query``.
        Must be nonnegative. A value of zero performs an exact search.

    p
        Minkowski distance norm passed to ``scipy.spatial.KDTree.query``.
        Must satisfy ``1 <= p <= infinity``. Common values are:

        - ``1``: Manhattan distance
        - ``2``: Euclidean distance
        - ``numpy.inf``: Chebyshev distance

    distance_upper_bound
        Maximum distance allowed for a returned neighbor. Must be positive
        or ``numpy.inf``.

    leafsize
        Number of observations at which KDTree switches to brute-force
        search. Must be a positive integer.

    Returns
    -------
    pandas.DataFrame
        Imputed matrix with the same index and columns as ``df``. Rows or
        columns that cannot be aligned with the reference are retained
        unchanged.

    Raises
    ------
    TypeError
        If ``df`` or ``ref`` is not a pandas DataFrame.

    ValueError
        If an argument is invalid, the reference contains missing values,
        labels are duplicated, or the matrices have no labels in common
        along the required dimension.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if not isinstance(ref, pd.DataFrame):
        raise TypeError("ref must be a pandas DataFrame.")

    if axis not in ("index", "columns"):
        raise ValueError("axis must be either 'index' or 'columns'.")

    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        raise ValueError("k must be a positive integer.")

    if not isinstance(eps, (int, float)) or isinstance(eps, bool) or eps < 0:
        raise ValueError("eps must be a nonnegative number.")

    if not isinstance(p, (int, float)) or isinstance(p, bool) or p < 1:
        raise ValueError("p must be greater than or equal to 1.")

    if (
        not isinstance(distance_upper_bound, (int, float))
        or isinstance(distance_upper_bound, bool)
        or distance_upper_bound <= 0
    ):
        raise ValueError(
            "distance_upper_bound must be positive or numpy.inf."
        )

    if (
        not isinstance(leafsize, int)
        or isinstance(leafsize, bool)
        or leafsize < 1
    ):
        raise ValueError("leafsize must be a positive integer.")

    if ref.isna().to_numpy().any():
        raise ValueError(
            "The external reference matrix must not contain missing values."
        )

    if df.index.has_duplicates:
        raise ValueError("The input matrix contains duplicated row labels.")

    if df.columns.has_duplicates:
        raise ValueError("The input matrix contains duplicated column labels.")

    if ref.index.has_duplicates:
        raise ValueError("The reference matrix contains duplicated row labels.")

    if ref.columns.has_duplicates:
        raise ValueError(
            "The reference matrix contains duplicated column labels."
        )

    result = df.copy()

    if axis == "columns":
        # Search reference columns using shared rows as features.
        common_labels = df.index[df.index.isin(ref.index)]

        if common_labels.empty:
            raise ValueError(
                "The input and reference matrices have no row labels in common."
            )

        input_aligned = df.loc[common_labels]
        ref_aligned = ref.loc[common_labels]

        n_candidates = ref_aligned.shape[1]
        if k > n_candidates:
            raise ValueError(
                f"k={k} exceeds the number of reference columns "
                f"available for neighbor searches ({n_candidates})."
            )

        # external_ref expects observations in rows. Transpose so that
        # columns/samples become observations.
        input_array = input_aligned.T.to_numpy(dtype=float)
        ref_array = ref_aligned.T.to_numpy(dtype=float)
        
        #fill NAs temporarily
        query_array, na_locations = _prepare_knn_query(
            input_array=input_array,
            ref_array=ref_array,
        )
        
        
        imputed_array = _external_ref(
            query_array,
            na_locations=na_locations,
            ref_data=ref_array,
            k=k,
            eps=eps,
            p=p,
            distance_upper_bound=distance_upper_bound,
            leafsize=leafsize,
        )

        imputed = pd.DataFrame(
            imputed_array.T,
            index=common_labels,
            columns=df.columns,
        )

        result.loc[common_labels, :] = imputed

    else:
        # Search reference rows using shared columns as features.
        common_labels = df.columns[df.columns.isin(ref.columns)]

        if common_labels.empty:
            raise ValueError(
                "The input and reference matrices have no column labels "
                "in common."
            )

        input_aligned = df.loc[:, common_labels]
        ref_aligned = ref.loc[:, common_labels]

        n_candidates = ref_aligned.shape[0]
        if k > n_candidates:
            raise ValueError(
                f"k={k} exceeds the number of reference rows "
                f"available for neighbor searches ({n_candidates})."
            )

        input_array = input_aligned.to_numpy(dtype=float)
        ref_array = ref_aligned.to_numpy(dtype=float)
        
        #fill NAs temporarily
        query_array, na_locations = _prepare_knn_query(
            input_array=input_array,
            ref_array=ref_array,
        )

        imputed_array = _external_ref(
            query_array,
            na_locations=na_locations,
            ref_data=ref_array,
            k=k,
            eps=eps,
            p=p,
            distance_upper_bound=distance_upper_bound,
            leafsize=leafsize,
        )

        imputed = pd.DataFrame(
            imputed_array,
            index=df.index,
            columns=common_labels,
        )

        result.loc[:, common_labels] = imputed

    return result


def toy_df(
    n_rows: int = 20,
    n_cols: int = 5,
    missingness: float = 0.2,
    min_val: float = 0.0,
    max_val: float = 1.0,
    rand_seed: int = 1234,
    sample_prefix: str | None = None,
) -> pd.DataFrame:
    """
    Generate a toy DataFrame containing randomly distributed missing values.

    Parameters
    ----------
    n_rows
        Number of rows.

    n_cols
        Number of columns.

    missingness
        Amount of missing data to introduce.

        - Values between 0 and 1 are interpreted as a fraction of all cells.
        - Values greater than or equal to 1 are interpreted as the number
          of missing cells.

    min_val
        Minimum generated value.

    max_val
        Maximum generated value.

    rand_seed
        Random seed.

    sample_prefix
        Prefix used for column names. If omitted, columns are assigned
        integer labels.

    Returns
    -------
    pandas.DataFrame
        Generated matrix containing exactly the requested number of
        missing values.
    """
    if not isinstance(n_rows, int) or isinstance(n_rows, bool) or n_rows < 1:
        raise ValueError("n_rows must be a positive integer.")

    if not isinstance(n_cols, int) or isinstance(n_cols, bool) or n_cols < 1:
        raise ValueError("n_cols must be a positive integer.")

    if missingness < 0:
        raise ValueError("missingness must be nonnegative.")

    if min_val >= max_val:
        raise ValueError("min_val must be smaller than max_val.")

    n_cells = n_rows * n_cols

    if missingness < 1:
        n_missing = int(missingness * n_cells)
    else:
        n_missing = int(missingness)

    if n_missing > n_cells:
        raise ValueError(
            f"Requested {n_missing} missing values, but the matrix "
            f"contains only {n_cells} cells."
        )

    rng = np.random.default_rng(rand_seed)

    values = rng.uniform(
        low=min_val,
        high=max_val,
        size=(n_rows, n_cols),
    )

    if n_missing > 0:
        missing_positions = rng.choice(
            n_cells,
            size=n_missing,
            replace=False,
        )
        values.flat[missing_positions] = np.nan

    if sample_prefix is None:
        columns = range(n_cols)
    else:
        columns = [
            f"{sample_prefix}_{i}"
            for i in range(n_cols)
        ]

    return pd.DataFrame(values, columns=columns)



Axis = Literal["index", "columns"]


def fill_moving_window(
    df: pd.DataFrame,
    axis: Axis = "index",
    naindex: int | None = None,
    wsize: int = 5,
    errors: Literal["raise", "coerce", "ignore"] = "coerce",
    func: Callable = np.mean) -> pd.DataFrame:
    """
    Impute missing values using a moving window.

    Missing values are replaced using a statistic calculated from
    neighboring values within a sliding window. The window is applied
    independently to each row or column.

    Parameters
    ----------
    df
        Input matrix with CpGs in rows and samples in columns.

    axis
        Dimension along which the moving window is applied.

        - ``"index"``: apply a horizontal moving window to each row (CpGs).
        - ``"columns"``: apply a vertical moving window to each column (samples).

    naindex
        Position of the missing value within the moving window.

        If ``None``, the missing value is placed at the center of the
        window whenever possible.

        Examples:

        - ``0``: use only values to the right.
        - ``-1``: use only values to the left.

    wsize
        Size of the moving window. The missing value itself is included
        when determining the window position.

    errors
        Behavior when the requested window extends beyond the available
        data.

        - ``"raise"``: raise an exception.
        - ``"coerce"``: automatically move the missing value toward the
          center of the window.
        - ``"ignore"``: leave the missing value unchanged.

    func
        Function used to summarize the values within each moving window.
        The default is ``numpy.mean``.

    Returns
    -------
    pandas.DataFrame
        Matrix with missing values imputed using moving-window statistics.

    Raises
    ------
    TypeError
        If ``df`` is not a pandas DataFrame.

    ValueError
        If an argument is invalid.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if axis not in ("index", "columns"):
        raise ValueError(
            "axis must be either 'index' or 'columns'."
        )

    if not isinstance(wsize, int) or isinstance(wsize, bool) or wsize < 1:
        raise ValueError("wsize must be a positive integer.")

    if errors not in ("raise", "coerce", "ignore"):
        raise ValueError(
            "errors must be one of 'raise', 'coerce', or 'ignore'."
        )

    if axis == "index":
        return _moving_window(
            df,
            nindex=naindex,
            wsize=wsize,
            errors=errors,
            func=SUMMARY_FUNCTIONS[func],
        )

    # axis == "columns"
    return _moving_window(
        df.T,
        nindex=naindex,
        wsize=wsize,
        errors=errors,
        func=SUMMARY_FUNCTIONS[func],
    ).T


def fill_knn(
    df: pd.DataFrame,
    axis: Axis = "columns",
    neighbors: int = 5,
    weights: Weights = "distance",
    keep_empty_features: bool = False,) -> pd.DataFrame:
    """
    Impute missing values using scikit-learn's KNNImputer.

    The input is assumed to be a DNA methylation matrix with CpGs in rows
    and samples in columns.

    Before KNN imputation:

    - Sample columns containing only missing values are always removed,
      because they contain no information for calculating distances.
    - CpG rows containing only missing values are removed by default.
      When ``keep_empty_features=True``, these CpGs are retained and filled
      with zeros after KNN imputation.

    Parameters
    ----------
    df
        Input DNA methylation matrix with CpGs in rows and samples in
        columns.

    axis
        Dimension in which nearest neighbors are searched.

        - ``"columns"``: search neighboring samples. This is the default.
        - ``"index"``: search neighboring CpGs.

    neighbors
        Number of nearest neighbors used for imputation. Must be a positive
        integer and cannot exceed the number of available observations along
        the selected search axis.

    weights
        Weighting method used to combine neighbor values.

        - ``"uniform"``: give every selected neighbor equal weight.
        - ``"distance"``: give closer neighbors greater weight.

    keep_empty_features
        How to handle CpG rows in which every sample value is missing.

        - ``False``: remove all-missing CpG rows.
        - ``True``: retain all-missing CpG rows and fill them with zeros.

        Sample columns containing only missing values are always removed,
        regardless of this setting.

    Returns
    -------
    pandas.DataFrame
        Imputed matrix. All-missing sample columns are removed. All-missing
        CpG rows are either removed or retained as zero-filled rows,
        depending on ``keep_empty_features``.

    Raises
    ------
    TypeError
        If an argument has an invalid type.

    ValueError
        If an argument has an invalid value or no usable data remain after
        removing empty samples and CpGs.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if df.empty:
        raise ValueError("df must not be empty.")

    if axis not in ("index", "columns"):
        raise ValueError("axis must be either 'index' or 'columns'.")

    if (
        not isinstance(neighbors, (int, np.integer))
        or isinstance(neighbors, bool)):
        raise TypeError("neighbors must be an integer.")

    if neighbors < 1:
        raise ValueError("neighbors must be a positive integer.")

    if weights not in ("uniform", "distance"):
        raise ValueError(
            "weights must be either 'uniform' or 'distance'.")

    if not isinstance(keep_empty_features, bool):
        raise TypeError("keep_empty_features must be a boolean.")

    try:
        numeric_df = df.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "df must contain only numeric values and missing values."
        ) from exc

    if np.isinf(numeric_df.to_numpy()).any():
        raise ValueError(
            "df contains infinite values. Replace or remove them before "
            "KNN imputation.")

    # Samples containing no observed CpGs cannot be used in KNN.
    empty_samples = numeric_df.isna().all(axis=0)
    retained_sample_names = numeric_df.columns[~empty_samples]

    working_df = numeric_df.loc[:, retained_sample_names]

    if working_df.shape[1] == 0:
        raise ValueError(
            "Every sample column contains only missing values.")

    # CpGs containing no observed values require explicit handling.
    empty_cpgs = working_df.isna().all(axis=1)
    empty_cpg_names = working_df.index[empty_cpgs]
    nonempty_cpg_names = working_df.index[~empty_cpgs]

    knn_df = working_df.loc[nonempty_cpg_names]

    if knn_df.shape[0] == 0:
        if keep_empty_features:
            return pd.DataFrame(
                0.0,
                index=working_df.index,
                columns=working_df.columns,)

        raise ValueError(
            "Every CpG row contains only missing values.")

    # KNNImputer always treats rows as observations and columns as features.
    if axis == "columns":
        # Samples are observations; CpGs are features.
        imputer_input = knn_df.T
    else:
        # CpGs are observations; samples are features.
        imputer_input = knn_df

    n_observations = imputer_input.shape[0]

    if neighbors > n_observations:
        raise ValueError(
            f"neighbors={neighbors} exceeds the number of observations "
            f"available for axis='{axis}' ({n_observations}).")

    imputer = KNNImputer(
        n_neighbors=neighbors,
        weights=weights,
        metric="nan_euclidean",)

    transformed = imputer.fit_transform(imputer_input)

    imputed_oriented = pd.DataFrame(
        transformed,
        index=imputer_input.index,
        columns=imputer_input.columns,)

    if axis == "columns":
        imputed_df = imputed_oriented.T
    else:
        imputed_df = imputed_oriented

    if keep_empty_features and len(empty_cpg_names) > 0:
        zero_rows = pd.DataFrame(
            0.0,
            index=empty_cpg_names,
            columns=working_df.columns,)

        imputed_df = pd.concat(
            [imputed_df, zero_rows],
            axis=0,)

        # Restore the original CpG order.
        imputed_df = imputed_df.loc[working_df.index]

    return imputed_df

def fill_buck(
    df: pd.DataFrame,
    eps: float = 1e-3,
    max_iter: int = 100,
    random_state: int | None = None,) -> pd.DataFrame:
    """
    Impute a DataFrame using iterative Buck regression imputation.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    imputed = _buck_iterative(
        data=df.to_numpy(dtype=float),
        eps=eps,
        max_iter=max_iter,
        random_state=random_state,)

    return pd.DataFrame(
        imputed,
        index=df.index,
        columns=df.columns,)

def fill_statistic(
    df: pd.DataFrame,
    axis: Axis = "index",
    statistic: Statistic = "mean",) -> pd.DataFrame:
    """
    Impute missing values using a summary statistic calculated from either
    rows (CpGs) or columns (samples).
    """
    if axis not in ("index", "columns"):
        raise ValueError("axis must be either 'index' or 'columns'.")

    try:
        func = STATISTICS[statistic]
    except KeyError:
        raise ValueError(
            "statistic must be one of 'mean', 'median', 'min', or 'max'."
        ) from None

    if axis == "columns":
        values = func(df)
        return df.fillna(values)

    # axis == "index"
    transposed = df.T
    values = func(transposed)
    return transposed.fillna(values).T


def fill_random(
    df: pd.DataFrame,
    axis: Axis = "index",) -> pd.DataFrame:
    """
    Impute missing values by randomly sampling observed values from the
    same row or column.

    Parameters
    ----------
    df
        Input matrix with CpGs in rows and samples in columns.

    axis
        Dimension from which replacement values are sampled.

        - ``"index"``: sample from the same row (CpG).
        - ``"columns"``: sample from the same column (sample).

    Returns
    -------
    pandas.DataFrame
        Imputed matrix.
    """
    if axis == "index":
        return _random_impute(df.T).T
    if axis == "columns":
        return _random_impute(df)
    raise ValueError("axis must be either 'index' or 'columns'.")

def fill_rf(
    df: pd.DataFrame,
    *,
    n_estimators: int = 100,
    max_iter: int = 10,
    eps: float = 1e-3,
    max_depth: int | None = None,
    min_samples_leaf: int = 1,
    random_state: int | None = 1234,
    n_jobs: int | None = -1,
    min_value: float = 0.0,
    max_value: float = 1.0,) -> pd.DataFrame:
    """
    Impute a DNA methylation matrix using iterative random forests.

    Rows are CpGs and columns are samples. Missing values are initialized
    using CpG row means. Each incomplete sample column is then predicted from
    the remaining sample columns across CpGs.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if df.empty:
        raise ValueError("df must not be empty.")

    try:
        values = df.to_numpy(dtype=float, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "df must contain only numeric and missing values."
        ) from exc

    imputed = _rf_iterative(
        values,
        n_estimators=n_estimators,
        max_iter=max_iter,
        eps=eps,
        max_depth=max_depth,
        min_samples_leaf=min_samples_leaf,
        random_state=random_state,
        n_jobs=n_jobs,
        min_value=min_value,
        max_value=max_value,)

    return pd.DataFrame(
        imputed,
        index=df.index.copy(),
        columns=df.columns.copy(),)


def fill_soft_impute(
    df: pd.DataFrame,
    *,
    shrinkage_value: float | None = None,
    convergence_threshold: float = 1e-3,
    max_iters: int = 100,
    max_rank: int | None = None,
    n_power_iterations: int = 1,
    init_fill_method: str = "zero",
    min_value: float | None = None,
    max_value: float | None = None,
    random_state: int | None = 1234,
    verbose: bool = False,) -> pd.DataFrame:
    """
    Impute missing values using SoftImpute matrix completion.

    SoftImpute repeatedly computes a low-rank singular-value decomposition,
    shrinks the singular values, reconstructs the matrix, and updates only
    the entries that were originally missing.

    The input is assumed to be a DNA methylation matrix with CpGs in rows and
    samples in columns.

    Parameters
    ----------
    df
        Input matrix with CpGs in rows and samples in columns.

    shrinkage_value
        Nonnegative amount subtracted from each singular value.

        If ``None``, the SoftImpute implementation selects a value from the
        initial matrix, typically based on its largest singular value.

    convergence_threshold
        Relative-change threshold used to stop iteration.

    max_iters
        Maximum number of SoftImpute iterations.

    max_rank
        Maximum number of singular components retained during each
        decomposition.

        If ``None``, the implementation may compute a full SVD. For large
        methylation matrices, a small value such as 5, 10, or 20 is strongly
        recommended.

    n_power_iterations
        Number of power iterations used by randomized SVD.

    init_fill_method
        Method used to initialize missing entries before SoftImpute begins.
        Supported values depend on the underlying ``Solver`` implementation.

    min_value
        Optional lower bound applied to imputed values. Use ``0.0`` for
        methylation beta values.

    max_value
        Optional upper bound applied to imputed values. Use ``1.0`` for
        methylation beta values.

    normalizer
        Optional normalizer accepted by the underlying SoftImpute class.

    random_state
        Random seed used by randomized SVD.

    verbose
        Whether the SoftImpute implementation should log progress.

    Returns
    -------
    pandas.DataFrame
        Imputed matrix with the original row and column labels preserved.

    Raises
    ------
    TypeError
        If ``df`` is not a pandas DataFrame.

    ValueError
        If the input matrix or parameters are invalid.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")
    if df.empty:
        raise ValueError("df must not be empty.")
    if df.index.has_duplicates:
        raise ValueError("df contains duplicated row labels.")
    if df.columns.has_duplicates:
        raise ValueError("df contains duplicated column labels.")
    if (
        not isinstance(convergence_threshold, (int, float, np.integer, np.floating))
        or isinstance(convergence_threshold, bool)
        or not np.isfinite(convergence_threshold)
        or convergence_threshold <= 0):
        raise ValueError(
            "convergence_threshold must be a finite number greater than zero.")

    if (
        not isinstance(max_iters, (int, np.integer))
        or isinstance(max_iters, bool)
        or max_iters < 1):
        raise ValueError("max_iters must be a positive integer.")

    if max_rank is not None and (
        not isinstance(max_rank, (int, np.integer))
        or isinstance(max_rank, bool)
        or max_rank < 1):
        raise ValueError("max_rank must be a positive integer or None.")

    if (
        not isinstance(n_power_iterations, (int, np.integer))
        or isinstance(n_power_iterations, bool)
        or n_power_iterations < 0):
        raise ValueError(
            "n_power_iterations must be a nonnegative integer."
        )

    if shrinkage_value is not None and (
        not isinstance(
            shrinkage_value,
            (int, float, np.integer, np.floating),)
        or isinstance(shrinkage_value, bool)
        or not np.isfinite(shrinkage_value)
        or shrinkage_value < 0):
        raise ValueError(
            "shrinkage_value must be a finite nonnegative number or None.")

    if min_value is not None and not np.isfinite(min_value):
        raise ValueError("min_value must be finite or None.")
    if max_value is not None and not np.isfinite(max_value):
        raise ValueError("max_value must be finite or None.")
    if (
        min_value is not None
        and max_value is not None
        and min_value >= max_value):
        raise ValueError("min_value must be smaller than max_value.")

    try:
        values = df.to_numpy(dtype=float, copy=True)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "df must contain only numeric values and missing values.") from exc

    if np.isinf(values).any():
        raise ValueError(
            "df contains infinite values. Replace or remove them before "
            "SoftImpute.")

    missing_mask = np.isnan(values)

    if not missing_mask.any():
        return df.astype(float).copy()

    all_missing_rows = missing_mask.all(axis=1)
    if all_missing_rows.any():
        rows = df.index[all_missing_rows].tolist()
        raise ValueError(
            "SoftImpute cannot reliably reconstruct rows containing only "
            f"missing values: {rows[:10]}"
            + (" ..." if len(rows) > 10 else ""))

    all_missing_columns = missing_mask.all(axis=0)
    if all_missing_columns.any():
        columns = df.columns[all_missing_columns].tolist()
        raise ValueError(
            "SoftImpute cannot reliably reconstruct columns containing only "
            f"missing values: {columns[:10]}"
            + (" ..." if len(columns) > 10 else ""))

    smaller_dimension = min(values.shape)

    if max_rank is not None and max_rank > smaller_dimension:
        raise ValueError(
            "max_rank cannot exceed the smaller matrix dimension "
            f"({smaller_dimension}); received {max_rank}.")

    imputer = SoftImpute(
        shrinkage_value=shrinkage_value,
        convergence_threshold=convergence_threshold,
        max_iters=max_iters,
        max_rank=max_rank,
        n_power_iterations=n_power_iterations,
        init_fill_method=init_fill_method,
        min_value=min_value,
        max_value=max_value,
        random_state=random_state,
        verbose=verbose,)

    imputed_values = imputer.fit_transform(values)

    if imputed_values.shape != values.shape:
        raise RuntimeError(
            "SoftImpute returned a matrix with an unexpected shape: "
            f"expected {values.shape}, received {imputed_values.shape}.")

    if not np.isfinite(imputed_values).all():
        raise RuntimeError(
            "SoftImpute produced missing or infinite values.")

    # Preserve observed values exactly in case the underlying implementation
    # introduces small numerical changes outside the missing positions.
    imputed_values[~missing_mask] = values[~missing_mask]

    return pd.DataFrame(
        imputed_values,
        index=df.index.copy(),
        columns=df.columns.copy(),)


def insert_na(
    df: pd.DataFrame,
    target_missing: int,
    *,
    random_state: int | np.random.Generator | None = 123,) -> pd.DataFrame:
    """
    Randomly insert missing values until a target total count is reached.

    Existing missing values are preserved. If the input already contains at
    least ``target_missing`` missing values, a copy of the original DataFrame
    is returned unchanged.

    Parameters
    ----------
    df
        Input DataFrame.

    target_missing
        Target total number of missing values in the returned DataFrame.

        For example, if ``target_missing=100`` and the input already contains
        20 missing values, 80 additional values are replaced with ``NaN``.

    random_state
        Integer seed, NumPy random generator, or ``None``. Use an integer for
        reproducible placement of missing values.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` containing the requested total number of missing values,
        unless the input already contains more missing values than requested.

    Raises
    ------
    TypeError
        If ``df`` is not a DataFrame, ``target_missing`` is not an integer,
        or ``random_state`` has an invalid type.

    ValueError
        If ``df`` is empty or ``target_missing`` is negative.

    RuntimeError
        If the final number of missing values is not the expected value.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {
    ...         "A": [0.0, 5.0, 10.0],
    ...         "B": [1.0, 6.0, 11.0],
    ...     }
    ... )
    >>> result = insert_na(df, target_missing=2, random_state=123)
    >>> int(result.isna().sum().sum())
    2
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if df.empty:
        raise ValueError("df must not be empty.")

    if (
        not isinstance(target_missing, Integral)
        or isinstance(target_missing, bool)
    ):
        raise TypeError("target_missing must be an integer.")

    if not (
        random_state is None
        or isinstance(random_state, Integral)
        or isinstance(random_state, np.random.Generator)
    ):
        raise TypeError(
            "random_state must be an integer, a NumPy Generator, or None."
        )

    if isinstance(random_state, bool):
        raise TypeError(
            "random_state must be an integer, a NumPy Generator, or None."
        )

    target_missing = int(target_missing)

    if target_missing < 0:
        raise ValueError("target_missing must be nonnegative.")

    result = df.copy(deep=True)

    missing_mask = result.isna().to_numpy()
    current_missing = int(missing_mask.sum())
    total_entries = result.size

    # Do not remove existing missing values.
    if target_missing <= current_missing:
        return result

    # A matrix cannot contain more missing values than total entries.
    target_missing = min(target_missing, total_entries)
    additional_needed = target_missing - current_missing

    if additional_needed == 0:
        return result

    observed_positions = np.flatnonzero(
        ~missing_mask.ravel())

    if additional_needed > observed_positions.size:
        raise ValueError(
            "There are not enough observed entries to reach the requested "
            f"target of {target_missing} missing values.")

    # Create or reuse the random number generator.
    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    selected_positions = rng.choice(
        observed_positions,
        size=additional_needed,
        replace=False,)

    values = result.to_numpy(copy=True)

    # Ensure the array can represent np.nan.
    if not (
        np.issubdtype(values.dtype, np.floating)
        or values.dtype == object):
        try:
            values = values.astype(float)
        except (TypeError, ValueError):
            values = values.astype(object)

    selected_rows, selected_columns = np.unravel_index(
        selected_positions,
        values.shape,)

    values[selected_rows, selected_columns] = np.nan

    result = pd.DataFrame(
        values,
        index=df.index.copy(),
        columns=df.columns.copy(),)

    final_missing = int(result.isna().to_numpy().sum())

    if final_missing != target_missing:
        raise RuntimeError(
            "Incorrect number of missing values after insertion: "
            f"expected {target_missing}, found {final_missing}.")

    return result


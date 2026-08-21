#!/usr/bin/env python3
#
# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: "deep.omics.lab@gmail.com"
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.


from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from numbers import Integral, Real
from typing import Literal

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import KNNImputer
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    r2_score,
)

LOGGER = logging.getLogger(__name__)

MorelModel = Literal["RF", "DNN"]


def _cluster_cols(data, k=2, random_state=0, n_init="auto"):
    """
    Binarize the columns of the DataFrame into two groups using K-means
    clustering (K=2). Note that the value of K is fixed and only supports
    two clusters.

    Missing values will be replaced with zero, while non-missing values will
    be replaced with one. K-means clustering (K=2) will then be applied to
    the columns to group samples with similar missing patterns.
    """
    #key is groupID, value is a list of samples (i.e., column names)
    group = {}
    df = data.copy()
    names = df.columns
    #replace non missing values with 1
    df = df.mask(df.notna(), 1)
    #replace missing values with 0
    df = df.fillna(0)
    #transpose df, since we want to cluster columns
    df = df.T

    # Initialize KMeans with desired number of clusters (k)
    kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=n_init).fit(df)
    labels = [str(i) for i in kmeans.labels_]

    for i,j in zip(names, labels):
        if j not in group:
            group[j] = [i]
        else:
            group[j].append(i)
    return group


def _validate_morel_groups(
    df: pd.DataFrame,
    groups: Mapping[str, Sequence[str]],
) -> tuple[list[str], list[str], list[str]]:
    """Validate a two-group sample definition."""
    if not isinstance(groups, Mapping):
        raise TypeError("group must be a mapping of group names to samples.")

    if len(groups) != 2:
        raise ValueError(
            "MOREL requires exactly two sample groups; "
            f"received {len(groups)}."
        )

    group_names = list(groups)

    group_0 = list(groups[group_names[0]])
    group_1 = list(groups[group_names[1]])

    if not group_0 or not group_1:
        raise ValueError("Each MOREL group must contain at least one sample.")

    all_samples = group_0 + group_1

    if len(all_samples) != len(set(all_samples)):
        duplicated = sorted(
            {
                sample
                for sample in all_samples
                if all_samples.count(sample) > 1
            }
        )
        raise ValueError(
            "Sample names must be unique across groups. "
            f"Duplicated samples: {duplicated}."
        )

    missing_samples = [
        sample for sample in all_samples
        if sample not in df.columns
    ]
    if missing_samples:
        raise ValueError(
            "The following grouped samples are absent from df: "
            f"{missing_samples}."
        )

    ungrouped_samples = [
        sample for sample in df.columns
        if sample not in all_samples
    ]
    if ungrouped_samples:
        raise ValueError(
            "Every DataFrame column must belong to one of the two groups. "
            f"Ungrouped samples: {ungrouped_samples}."
        )

    return group_names, group_0, group_1


def _knn_impute_sporadic(
    df: pd.DataFrame,
    *,
    neighbors: int,
    weights: Literal["uniform", "distance"],
) -> pd.DataFrame:
    """
    Impute sporadic missing values while treating samples as observations.

    The input has CpGs in rows and samples in columns. It is transposed before
    KNN imputation so that samples are observations and CpGs are features.
    """
    if not df.isna().to_numpy().any():
        return df.copy()

    if df.isna().all(axis=1).any():
        rows = df.index[df.isna().all(axis=1)].tolist()
        raise ValueError(
            "KNN preprocessing cannot impute CpGs that are entirely missing "
            f"in the supplied matrix: {rows[:10]}"
            + (" ..." if len(rows) > 10 else "")
        )

    if df.isna().all(axis=0).any():
        columns = df.columns[df.isna().all(axis=0)].tolist()
        raise ValueError(
            "KNN preprocessing cannot use samples that are entirely missing: "
            f"{columns}."
        )

    n_samples = df.shape[1]
    effective_neighbors = min(neighbors, n_samples)

    if effective_neighbors < 1:
        raise ValueError(
            "At least one sample is required for KNN preprocessing."
        )

    imputer = KNNImputer(
        n_neighbors=effective_neighbors,
        weights=weights,
        keep_empty_features=False,
    )

    transformed = imputer.fit_transform(df.T)

    if transformed.shape[1] != df.shape[0]:
        raise RuntimeError(
            "KNN preprocessing changed the number of CpG features."
        )

    return pd.DataFrame(
        transformed.T,
        index=df.index.copy(),
        columns=df.columns.copy(),
    )


def _predict_block_rf(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_predict: pd.DataFrame,
    *,
    n_iter: int,
    train_size: float,
    random_state: int,
    n_estimators: int,
    max_depth: int | None,
    min_samples_leaf: int,
    max_features: float | str | None,
    n_jobs: int | None,
    group_name: str,
) -> np.ndarray:
    """Predict a block of missing values using repeated random forests."""
    if x_predict.empty:
        return np.empty((0, y_train.shape[1]), dtype=float)

    predictions: list[np.ndarray] = []

    for iteration in range(n_iter):
        split_seed = random_state + iteration

        x_fit, x_test, y_fit, y_test = train_test_split(
            x_train,
            y_train,
            train_size=train_size,
            random_state=split_seed,
        )

        model = RandomForestRegressor(
            n_estimators=n_estimators,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            max_features=max_features,
            random_state=split_seed,
            n_jobs=n_jobs,
        )

        model.fit(x_fit, y_fit)

        test_prediction = model.predict(x_test)

        score = r2_score(
            y_test,
            test_prediction,
            multioutput="uniform_average",
        )

        LOGGER.info(
            'RF iteration %d/%d for group "%s": test R2 = %.4f.',
            iteration + 1,
            n_iter,
            group_name,
            score,
        )

        # Refit using all complete training rows before predicting the
        # block-missing rows.
        model.fit(x_train, y_train)
        predictions.append(model.predict(x_predict))

    return np.mean(predictions, axis=0)


def _predict_block_dnn(
    x_train: pd.DataFrame,
    y_train: pd.DataFrame,
    x_predict: pd.DataFrame,
    *,
    epochs: int,
    train_size: float,
    random_state: int,
    group_name: str,
) -> np.ndarray:
    """Predict a block of missing values using a dense neural network."""
    if x_predict.empty:
        return np.empty((0, y_train.shape[1]), dtype=float)

    try:
        import tensorflow as tf
        from tensorflow.keras import Sequential
        from tensorflow.keras.callbacks import EarlyStopping
        from tensorflow.keras.layers import Dense, Input
    except ImportError as exc:
        raise ImportError(
            "TensorFlow is required when model='DNN'. Install it with "
            "'pip install tensorflow'."
        ) from exc

    np.random.seed(random_state)
    tf.keras.utils.set_random_seed(random_state)

    x_fit, x_test, y_fit, y_test = train_test_split(
        x_train.to_numpy(dtype=float),
        y_train.to_numpy(dtype=float),
        train_size=train_size,
        random_state=random_state,
    )

    input_dimension = x_train.shape[1]
    output_dimension = y_train.shape[1]
    hidden_units = max(
        2,
        int((input_dimension + output_dimension) / 2) + 1,
    )

    def build_model() -> Sequential:
        network = Sequential(
            [
                Input(shape=(input_dimension,)),
                Dense(
                    hidden_units,
                    kernel_initializer="he_uniform",
                    activation="relu",
                ),
                Dense(
                    output_dimension,
                    kernel_initializer="he_uniform",
                    activation="linear",
                ),
            ]
        )

        network.compile(
            optimizer="adam",
            loss="mae",
            metrics=["root_mean_squared_error"],
        )
        return network

    model = build_model()

    early_stopping = EarlyStopping(
        monitor="val_loss",
        patience=10,
        restore_best_weights=True,
    )

    model.fit(
        x_fit,
        y_fit,
        validation_data=(x_test, y_test),
        epochs=epochs,
        verbose=0,
        callbacks=[early_stopping],
    )

    test_prediction = model.predict(x_test, verbose=0)

    mae = mean_absolute_error(y_test, test_prediction)
    rmse = np.sqrt(mean_squared_error(y_test, test_prediction))
    score = r2_score(
        y_test,
        test_prediction,
        multioutput="uniform_average",
    )

    LOGGER.info(
        'DNN validation for group "%s": MAE=%.6f, RMSE=%.6f, R2=%.4f.',
        group_name,
        mae,
        rmse,
        score,
    )

    # Fit a new model on all available complete rows.
    tf.keras.backend.clear_session()
    tf.keras.utils.set_random_seed(random_state)

    final_model = build_model()
    final_model.fit(
        x_train.to_numpy(dtype=float),
        y_train.to_numpy(dtype=float),
        epochs=epochs,
        verbose=0,
    )

    return final_model.predict(
        x_predict.to_numpy(dtype=float),
        verbose=0,
    )


def fill_morel(
    df: pd.DataFrame,
    group: Mapping[str, Sequence[str]] | None = None,
    *,
    model: MorelModel = "RF",
    decimal: int = 5,
    knn_neighbors: int = 5,
    knn_weights: Literal["uniform", "distance"] = "uniform",
    n_iter: int = 10,
    random_state: int = 100,
    n_jobs: int | None = -1,
    train_size: float = 0.75,
    n_estimators: int = 100,
    max_depth: int | None = 30,
    min_samples_leaf: int = 1,
    max_features: float | str | None = 1.0,
    epochs: int = 100,
    min_value: float | None = 0.0,
    max_value: float | None = 1.0,
) -> pd.DataFrame:
    """
    Impute block-wise missing values using the MOREL strategy.

    The input is expected to be a DNA methylation matrix with CpGs in rows
    and samples in columns. Samples are divided into exactly two groups.

    CpGs that are complete in both groups are used to train a model relating
    one sample group to the other. When a CpG is entirely missing in one
    group but observed in the other, its missing group is predicted from its
    observed group.

    Sporadic missing values are first imputed using K-nearest neighbors.
    Block-wise missing values are then predicted using Random Forest or a
    deep neural network.

    Parameters
    ----------
    df
        DNA methylation matrix with CpGs in rows and samples in columns.

    group
        Mapping containing exactly two sample groups. For example::

            {
                "A": ["sample1", "sample2", "sample3"],
                "B": ["sample4", "sample5", "sample6"],
            }

        Every column in ``df`` must occur in exactly one group. If ``None``,
        ``cluster_cols(df)`` is used to identify two groups automatically.

    model
        Secondary model used to predict block-wise missing values:

        - ``"RF"``: Random Forest regression.
        - ``"DNN"``: dense neural network.

    decimal
        Number of decimal places used to round the final matrix.

    knn_neighbors
        Number of neighbors used to impute sporadic missing values.

    knn_weights
        KNN weighting strategy: ``"uniform"`` or ``"distance"``.

    n_iter
        Number of repeated train/test splits and forest fits when
        ``model="RF"``. Predictions from all repetitions are averaged.

    random_state
        Random seed used for data splitting and model fitting.

    n_jobs
        Number of parallel jobs used by Random Forest. ``-1`` uses all
        available processors.

    train_size
        Fraction of complete CpGs used for model training during validation.

    n_estimators
        Number of trees in each Random Forest.

    max_depth
        Maximum depth of each Random Forest tree.

    min_samples_leaf
        Minimum number of observations required in a terminal forest leaf.

    max_features
        Number or fraction of predictors considered at each forest split.

    epochs
        Maximum number of training epochs when ``model="DNN"``.

    min_value
        Optional lower bound applied to imputed values. Use ``0.0`` for
        methylation beta values.

    max_value
        Optional upper bound applied to imputed values. Use ``1.0`` for
        methylation beta values.

    Returns
    -------
    pandas.DataFrame
        Imputed matrix with the original index and column order preserved.

    Raises
    ------
    TypeError
        If an input has an invalid type.

    ValueError
        If a parameter, group definition, or input matrix is invalid.

    Notes
    -----
    CpGs that are missing across every sample cannot be predicted and remain
    missing in the returned DataFrame.

    This implementation is intended for systematic block-wise missingness.
    Users may alternatively impute sporadic missing values before calling
    this function.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if df.empty:
        raise ValueError("df must not be empty.")

    if df.index.has_duplicates:
        raise ValueError("df contains duplicated CpG identifiers.")

    if df.columns.has_duplicates:
        raise ValueError("df contains duplicated sample names.")

    if model not in ("RF", "DNN"):
        raise ValueError("model must be either 'RF' or 'DNN'.")

    if (
        not isinstance(decimal, Integral)
        or isinstance(decimal, bool)
        or decimal < 0
    ):
        raise ValueError("decimal must be a nonnegative integer.")

    if (
        not isinstance(knn_neighbors, Integral)
        or isinstance(knn_neighbors, bool)
        or knn_neighbors < 1
    ):
        raise ValueError("knn_neighbors must be a positive integer.")

    if knn_weights not in ("uniform", "distance"):
        raise ValueError(
            "knn_weights must be either 'uniform' or 'distance'."
        )

    if (
        not isinstance(n_iter, Integral)
        or isinstance(n_iter, bool)
        or n_iter < 1
    ):
        raise ValueError("n_iter must be a positive integer.")

    if (
        not isinstance(train_size, Real)
        or isinstance(train_size, bool)
        or not 0.0 < train_size < 1.0
    ):
        raise ValueError("train_size must be between 0 and 1.")

    if (
        not isinstance(epochs, Integral)
        or isinstance(epochs, bool)
        or epochs < 1
    ):
        raise ValueError("epochs must be a positive integer.")

    if (
        min_value is not None
        and max_value is not None
        and min_value >= max_value
    ):
        raise ValueError("min_value must be smaller than max_value.")

    try:
        values = df.astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "df must contain only numeric values and missing values."
        ) from exc

    if np.isinf(values.to_numpy()).any():
        raise ValueError("df must not contain infinite values.")

    if group is None:
        LOGGER.info(
            "No groups supplied; identifying two sample groups using K-means."
        )
        group = _cluster_cols(values)

    group_names, group_0, group_1 = _validate_morel_groups(
        values,
        group,
    )

    LOGGER.info(
        'Group "%s" contains %d samples.',
        group_names[0],
        len(group_0),
    )
    LOGGER.info(
        'Group "%s" contains %d samples.',
        group_names[1],
        len(group_1),
    )

    # CpGs missing across every sample cannot be modeled.
    globally_missing = values.isna().all(axis=1)
    globally_missing_df = values.loc[globally_missing].copy()

    usable = values.loc[~globally_missing].copy()

    group_0_all_missing = usable[group_0].isna().all(axis=1)
    group_1_all_missing = usable[group_1].isna().all(axis=1)

    LOGGER.info(
        '%d CpGs are entirely missing in group "%s".',
        int(group_0_all_missing.sum()),
        group_names[0],
    )
    LOGGER.info(
        '%d CpGs are entirely missing in group "%s".',
        int(group_1_all_missing.sum()),
        group_names[1],
    )

    # A non-global row cannot be completely missing in both groups.
    overlap = group_0_all_missing & group_1_all_missing
    if overlap.any():
        raise RuntimeError(
            "One or more non-global rows were classified as entirely "
            "missing in both groups."
        )

    # Rows suitable for training have at least one observed value in each
    # group. Any sporadic missing entries are handled using KNN.
    training_rows = ~(group_0_all_missing | group_1_all_missing)
    training_df = usable.loc[training_rows].copy()

    sporadic_count = int(training_df.isna().to_numpy().sum())
    if sporadic_count:
        LOGGER.info(
            "Imputing %d sporadic missing values with KNN.",
            sporadic_count,
        )
        training_df = _knn_impute_sporadic(
            training_df,
            neighbors=knn_neighbors,
            weights=knn_weights,
        )

    group_0_train = training_df[group_0]
    group_1_train = training_df[group_1]

    if len(training_df) < 3:
        raise ValueError(
            "At least three complete or KNN-completable CpG rows are "
            "required to train the secondary model."
        )

    # Rows where group 0 is the missing target and group 1 is the predictor.
    group_0_block = usable.loc[group_0_all_missing].copy()

    if not group_0_block.empty:
        predictor_1 = group_0_block[group_1]

        predictor_missing = int(predictor_1.isna().to_numpy().sum())
        if predictor_missing:
            LOGGER.info(
                "Imputing %d sporadic predictor values in group \"%s\" "
                "with KNN.",
                predictor_missing,
                group_names[1],
            )
            predictor_1 = _knn_impute_sporadic(
                predictor_1,
                neighbors=knn_neighbors,
                weights=knn_weights,
            )
            group_0_block.loc[:, group_1] = predictor_1

        LOGGER.info(
            'Predicting %d block-missing CpGs in group "%s" using %s.',
            len(group_0_block),
            group_names[0],
            model,
        )

        if model == "RF":
            prediction = _predict_block_rf(
                group_1_train,
                group_0_train,
                predictor_1,
                n_iter=n_iter,
                train_size=train_size,
                random_state=random_state,
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_leaf=min_samples_leaf,
                max_features=max_features,
                n_jobs=n_jobs,
                group_name=group_names[0],
            )
        else:
            prediction = _predict_block_dnn(
                group_1_train,
                group_0_train,
                predictor_1,
                epochs=epochs,
                train_size=train_size,
                random_state=random_state,
                group_name=group_names[0],
            )

        group_0_block.loc[:, group_0] = prediction

    # Rows where group 1 is the missing target and group 0 is the predictor.
    group_1_block = usable.loc[group_1_all_missing].copy()

    if not group_1_block.empty:
        predictor_0 = group_1_block[group_0]

        predictor_missing = int(predictor_0.isna().to_numpy().sum())
        if predictor_missing:
            LOGGER.info(
                "Imputing %d sporadic predictor values in group \"%s\" "
                "with KNN.",
                predictor_missing,
                group_names[0],
            )
            predictor_0 = _knn_impute_sporadic(
                predictor_0,
                neighbors=knn_neighbors,
                weights=knn_weights,
            )
            group_1_block.loc[:, group_0] = predictor_0

        LOGGER.info(
            'Predicting %d block-missing CpGs in group "%s" using %s.',
            len(group_1_block),
            group_names[1],
            model,
        )

        if model == "RF":
            prediction = _predict_block_rf(
                group_0_train,
                group_1_train,
                predictor_0,
                n_iter=n_iter,
                train_size=train_size,
                random_state=random_state,
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_leaf=min_samples_leaf,
                max_features=max_features,
                n_jobs=n_jobs,
                group_name=group_names[1],
            )
        else:
            prediction = _predict_block_dnn(
                group_0_train,
                group_1_train,
                predictor_0,
                epochs=epochs,
                train_size=train_size,
                random_state=random_state,
                group_name=group_names[1],
            )

        group_1_block.loc[:, group_1] = prediction

    result = pd.concat(
        [
            training_df,
            group_0_block,
            group_1_block,
            globally_missing_df,
        ],
        axis=0,
    )

    result = result.reindex(
        index=df.index,
        columns=df.columns,
    )

    if min_value is not None or max_value is not None:
        result = result.clip(
            lower=min_value,
            upper=max_value,
        )

    return result.round(decimal)

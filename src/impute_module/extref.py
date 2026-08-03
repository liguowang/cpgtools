from typing import Literal

import numpy as np
import pandas as pd


SearchAxis = Literal["index", "columns"]


def fill_ref(
    self,
    ref: pd.DataFrame,
    axis: SearchAxis = "columns",
    k: int = 3,
    eps: float = 0.0,
    p: float = 2.0,
    distance_upper_bound: float = np.inf,
    leafsize: int = 10,
) -> pd.DataFrame:
    """
    Impute missing values using nearest neighbors from an external reference.

    This method is similar to fast KNN imputation, except that nearest
    neighbors are searched in an external reference matrix rather than in
    the input matrix. Missing values are replaced using the corresponding
    values from the selected reference neighbors.

    Parameters
    ----------
    ref
        Complete external reference matrix. It must not contain missing
        values.

        When ``axis="columns"``, columns are treated as observations and
        rows as features. The input and reference matrices are aligned by
        their row labels.

        When ``axis="index"``, rows are treated as observations and columns
        as features. The input and reference matrices are aligned by their
        column labels.

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
        Imputed matrix with the same index and columns as the input matrix.
        Rows or columns that cannot be aligned with the reference are
        retained unchanged.

    Raises
    ------
    TypeError
        If ``ref`` is not a pandas DataFrame.

    ValueError
        If an argument is invalid, the reference contains missing values,
        labels are duplicated, or the input and reference matrices have no
        labels in common along the required dimension.
    """
    if not isinstance(ref, pd.DataFrame):
        raise TypeError("ref must be a pandas DataFrame.")

    if axis not in ("index", "columns"):
        raise ValueError("axis must be either 'index' or 'columns'.")

    if not isinstance(k, int) or isinstance(k, bool) or k < 1:
        raise ValueError("k must be a positive integer.")

    if eps < 0:
        raise ValueError("eps must be nonnegative.")

    if p < 1:
        raise ValueError("p must be greater than or equal to 1.")

    if distance_upper_bound <= 0:
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

    if self.df.index.has_duplicates:
        raise ValueError("The input matrix contains duplicated row labels.")

    if self.df.columns.has_duplicates:
        raise ValueError("The input matrix contains duplicated column labels.")

    if ref.index.has_duplicates:
        raise ValueError("The reference matrix contains duplicated row labels.")

    if ref.columns.has_duplicates:
        raise ValueError(
            "The reference matrix contains duplicated column labels."
        )

    result = self.df.copy()

    if axis == "columns":
        # Search for neighboring reference columns using shared rows as
        # features. Preserve the input row order.
        common_labels = self.df.index[
            self.df.index.isin(ref.index)
        ]

        if common_labels.empty:
            raise ValueError(
                "The input and reference matrices have no row labels "
                "in common."
            )

        input_aligned = self.df.loc[common_labels]
        ref_aligned = ref.loc[common_labels]

        n_candidates = ref_aligned.shape[1]
        if k > n_candidates:
            raise ValueError(
                f"k={k} exceeds the number of reference columns "
                f"available for neighbor searches ({n_candidates})."
            )

        # external_ref expects observations in rows, so transpose both
        # matrices before searching columns.
        input_array = input_aligned.T.to_numpy()
        ref_array = ref_aligned.T.to_numpy()
        na_locations = nan_indices(input_array)

        imputed_array = external_ref(
            input_array,
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
            columns=self.df.columns,
        )

        result.loc[common_labels, :] = imputed

    else:
        # Search for neighboring reference rows using shared columns as
        # features. Preserve the input column order.
        common_labels = self.df.columns[
            self.df.columns.isin(ref.columns)
        ]

        if common_labels.empty:
            raise ValueError(
                "The input and reference matrices have no column labels "
                "in common."
            )

        input_aligned = self.df.loc[:, common_labels]
        ref_aligned = ref.loc[:, common_labels]

        n_candidates = ref_aligned.shape[0]
        if k > n_candidates:
            raise ValueError(
                f"k={k} exceeds the number of reference rows "
                f"available for neighbor searches ({n_candidates})."
            )

        input_array = input_aligned.to_numpy()
        ref_array = ref_aligned.to_numpy()
        na_locations = nan_indices(input_array)

        imputed_array = external_ref(
            input_array,
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
            index=self.df.index,
            columns=common_labels,
        )

        result.loc[:, common_labels] = imputed

    return result

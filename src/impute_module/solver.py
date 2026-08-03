# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import annotations

import warnings
from abc import ABC, abstractmethod
from collections.abc import Callable
from numbers import Integral, Real
from typing import Any, Literal
import numpy as np
from numpy.typing import ArrayLike, NDArray
from sklearn.utils import check_array
from .common import generate_random_column_samples

FillMethod = Literal["zero", "mean", "median", "min", "random"]
FloatArray = NDArray[np.float64]
BoolArray = NDArray[np.bool_]

class Solver(ABC):
    """
    Base class for matrix-completion and missing-value imputation algorithms.

    Subclasses must implement :meth:`solve`, which receives an initialized
    numeric matrix and a Boolean mask identifying entries that were originally
    missing.
    """

    VALID_FILL_METHODS: tuple[FillMethod, ...] = (
        "zero",
        "mean",
        "median",
        "min",
        "random",
    )

    def __init__(
        self,
        fill_method: FillMethod = "zero",
        min_value: float | None = None,
        max_value: float | None = None,
        normalizer: Any | None = None,
        random_state: int | None = 1234,
    ) -> None:
        """
        Parameters
        ----------
        fill_method
            Method used to initialize missing entries before calling
            :meth:`solve`.

        min_value
            Optional lower bound applied to the completed matrix.

        max_value
            Optional upper bound applied to the completed matrix.

        normalizer
            Optional object providing ``fit_transform`` and
            ``inverse_transform`` methods.

        random_state
            Seed used by random initialization. Use ``None`` for
            nondeterministic sampling.
        """
        self._validate_constructor_parameters(
            fill_method=fill_method,
            min_value=min_value,
            max_value=max_value,
            normalizer=normalizer,
            random_state=random_state,
        )

        self.fill_method = fill_method
        self.min_value = (
            None if min_value is None else float(min_value)
        )
        self.max_value = (
            None if max_value is None else float(max_value)
        )
        self.normalizer = normalizer
        self.random_state = random_state

    @classmethod
    def _validate_constructor_parameters(
        cls,
        *,
        fill_method: str,
        min_value: float | None,
        max_value: float | None,
        normalizer: Any | None,
        random_state: int | None,
    ) -> None:
        """Validate parameters shared by all solver subclasses."""
        if fill_method not in cls.VALID_FILL_METHODS:
            choices = ", ".join(repr(x) for x in cls.VALID_FILL_METHODS)
            raise ValueError(
                f"fill_method must be one of {choices}; "
                f"received {fill_method!r}."
            )

        if min_value is not None and (
            not isinstance(min_value, Real)
            or isinstance(min_value, bool)
            or not np.isfinite(min_value)
        ):
            raise ValueError(
                "min_value must be a finite number or None."
            )

        if max_value is not None and (
            not isinstance(max_value, Real)
            or isinstance(max_value, bool)
            or not np.isfinite(max_value)
        ):
            raise ValueError(
                "max_value must be a finite number or None."
            )

        if (
            min_value is not None
            and max_value is not None
            and min_value >= max_value
        ):
            raise ValueError(
                "min_value must be smaller than max_value."
            )

        if random_state is not None and (
            not isinstance(random_state, Integral)
            or isinstance(random_state, bool)
        ):
            raise TypeError(
                "random_state must be an integer or None."
            )

        if normalizer is not None:
            required_methods = (
                "fit_transform",
                "inverse_transform",
            )

            missing_methods = [
                method
                for method in required_methods
                if not callable(getattr(normalizer, method, None))
            ]

            if missing_methods:
                raise TypeError(
                    "normalizer must provide callable "
                    f"{', '.join(missing_methods)} method(s)."
                )

    def __repr__(self) -> str:
        """Return an unambiguous representation of the solver."""
        fields: list[str] = []

        for key, value in sorted(self.__dict__.items()):
            if value is None or isinstance(
                value,
                (str, int, float, bool),
            ):
                fields.append(f"{key}={value!r}")

        return (
            f"{self.__class__.__name__}"
            f"({', '.join(fields)})"
        )

    __str__ = __repr__

    @staticmethod
    def _validate_missing_mask(
        missing_mask: ArrayLike,
        expected_shape: tuple[int, int],
    ) -> BoolArray:
        """Validate and return a Boolean missing-value mask."""
        mask = np.asarray(missing_mask)

        if mask.ndim != 2:
            raise ValueError(
                "missing_mask must be two-dimensional."
            )

        if mask.shape != expected_shape:
            raise ValueError(
                "missing_mask must have the same shape as the matrix; "
                f"expected {expected_shape}, received {mask.shape}."
            )

        if mask.dtype != np.bool_:
            raise TypeError(
                "missing_mask must contain Boolean values."
            )

        return mask

    @staticmethod
    def _check_missing_value_mask(
        missing_mask: BoolArray,
    ) -> None:
        """Validate that the matrix is suitable for imputation."""
        if not missing_mask.any():
            warnings.warn(
                "The input matrix does not contain missing values.",
                UserWarning,
                stacklevel=3,
            )

        if missing_mask.all():
            raise ValueError(
                "The input matrix contains no observed values."
            )

    @staticmethod
    def _validate_fill_values(
        fill_values: ArrayLike | float,
        expected_count: int,
        column_index: int,
    ) -> FloatArray:
        """Normalize and validate values returned by a fill function."""
        values = np.asarray(fill_values, dtype=float)

        if values.ndim == 0:
            values = np.full(
                expected_count,
                float(values),
                dtype=float,
            )
        else:
            values = values.reshape(-1)

        if values.size != expected_count:
            raise ValueError(
                f"Fill function returned {values.size} value(s) for "
                f"column {column_index}, but {expected_count} were "
                "required."
            )

        if not np.isfinite(values).all():
            raise ValueError(
                f"Could not generate finite initialization values for "
                f"column {column_index}."
            )

        return values

    def _fill_columns_with_fn(
        self,
        values: FloatArray,
        missing_mask: BoolArray,
        column_function: Callable[[FloatArray], ArrayLike | float],
    ) -> None:
        """Fill missing values independently within each column."""
        for column_index in range(values.shape[1]):
            column_missing = missing_mask[:, column_index]
            number_missing = int(column_missing.sum())

            if number_missing == 0:
                continue

            column = values[:, column_index]
            observed = column[~column_missing]

            if observed.size == 0:
                raise ValueError(
                    "Cannot initialize column "
                    f"{column_index} because it contains only "
                    "missing values."
                )

            fill_values = column_function(column)

            values[column_missing, column_index] = (
                self._validate_fill_values(
                    fill_values,
                    expected_count=number_missing,
                    column_index=column_index,
                )
            )

    def fill(
        self,
        X: ArrayLike,
        missing_mask: ArrayLike,
        fill_method: FillMethod | None = None,
        *,
        inplace: bool = False,
    ) -> FloatArray:
        """
        Initialize missing entries in a numeric matrix.

        Parameters
        ----------
        X
            Two-dimensional numeric matrix.

        missing_mask
            Boolean matrix identifying entries that were originally missing.

        fill_method
            Initialization strategy. When omitted, ``self.fill_method`` is
            used.

            ``"zero"``
                Fill missing entries with zero.

            ``"mean"``
                Fill using the corresponding column mean.

            ``"median"``
                Fill using the corresponding column median.

            ``"min"``
                Fill using the corresponding column minimum.

            ``"random"``
                Draw values from a normal distribution estimated from the
                corresponding column.

        inplace
            Modify the supplied floating-point NumPy array when possible.
            Otherwise, return a copy.

        Returns
        -------
        numpy.ndarray
            Initialized floating-point matrix.
        """
        checked = check_array(
            X,
            dtype=np.float64,
            ensure_2d=True,
            copy=not inplace,
            ensure_all_finite=False,
        )

        mask = self._validate_missing_mask(
            missing_mask,
            checked.shape,
        )

        method = (
            self.fill_method
            if fill_method is None
            else fill_method
        )

        if method not in self.VALID_FILL_METHODS:
            choices = ", ".join(
                repr(x) for x in self.VALID_FILL_METHODS
            )
            raise ValueError(
                f"fill_method must be one of {choices}; "
                f"received {method!r}."
            )

        values = checked if inplace else checked.copy()

        if not mask.any():
            return values

        if method == "zero":
            values[mask] = 0.0

        elif method == "mean":
            self._fill_columns_with_fn(
                values,
                mask,
                np.nanmean,
            )

        elif method == "median":
            self._fill_columns_with_fn(
                values,
                mask,
                np.nanmedian,
            )

        elif method == "min":
            self._fill_columns_with_fn(
                values,
                mask,
                np.nanmin,
            )

        elif method == "random":
            rng = np.random.default_rng(self.random_state)

            self._fill_columns_with_fn(
                values,
                mask,
                lambda column: generate_random_column_samples(
                    column,
                    random_state=rng,
                ),
            )

        if not np.isfinite(values).all():
            raise RuntimeError(
                "Initial filling produced missing or infinite values."
            )

        return values

    def prepare_input_data(
        self,
        X: ArrayLike,
    ) -> tuple[FloatArray, BoolArray]:
        """
        Validate input and return a floating-point matrix and missing mask.
        """
        values = check_array(
            X,
            dtype=np.float64,
            ensure_2d=True,
            copy=True,
            ensure_all_finite=False,
        )

        if np.isinf(values).any():
            raise ValueError(
                "The input matrix contains infinite values."
            )

        missing_mask = np.isnan(values)
        self._check_missing_value_mask(missing_mask)

        return values, missing_mask

    def clip(
        self,
        X: ArrayLike,
        *,
        inplace: bool = False,
    ) -> FloatArray:
        """
        Clip values to the configured global lower and upper bounds.
        """
        values = np.asarray(X, dtype=float)

        if not inplace:
            values = values.copy()

        return np.clip(
            values,
            a_min=self.min_value,
            a_max=self.max_value,
            out=values,
        )

    def project_result(
        self,
        X: ArrayLike,
    ) -> FloatArray:
        """
        Undo optional normalization and apply configured value bounds.
        """
        values = np.asarray(X, dtype=float)

        if values.ndim != 2:
            raise ValueError(
                "The completed matrix must be two-dimensional."
            )

        if self.normalizer is not None:
            values = np.asarray(
                self.normalizer.inverse_transform(values),
                dtype=float,
            )

        values = self.clip(values)

        if not np.isfinite(values).all():
            raise RuntimeError(
                "Projection produced missing or infinite values."
            )

        return values

    @abstractmethod
    def solve(
        self,
        X: FloatArray,
        missing_mask: BoolArray,
    ) -> FloatArray:
        """
        Complete an initialized matrix.

        Subclasses must return a NumPy array with the same shape as ``X``.
        """
        raise NotImplementedError

    def fit_transform(
        self,
        X: ArrayLike,
        y: Any | None = None,
    ) -> FloatArray:
        """
        Initialize, solve, reverse normalization, and restore observed values.
        """
        del y  # This base interface does not use supervised targets.

        original, missing_mask = self.prepare_input_data(X)

        if not missing_mask.any():
            return original.copy()

        observed_mask = ~missing_mask
        working = original.copy()

        if self.normalizer is not None:
            working = np.asarray(
                self.normalizer.fit_transform(working),
                dtype=float,
            )

            if working.shape != original.shape:
                raise ValueError(
                    "normalizer.fit_transform() changed the matrix shape: "
                    f"expected {original.shape}, received {working.shape}."
                )

        initialized = self.fill(
            working,
            missing_mask,
            inplace=True,
        )

        completed = self.solve(
            initialized,
            missing_mask,
        )

        if not isinstance(completed, np.ndarray):
            raise TypeError(
                f"{self.__class__.__name__}.solve() must return a NumPy "
                f"array, not {type(completed).__name__}."
            )

        if completed.shape != original.shape:
            raise ValueError(
                f"{self.__class__.__name__}.solve() returned shape "
                f"{completed.shape}; expected {original.shape}."
            )

        result = self.project_result(completed)

        # Preserve the original observed values exactly.
        result[observed_mask] = original[observed_mask]

        if np.isnan(result[missing_mask]).any():
            raise RuntimeError(
                f"{self.__class__.__name__} left one or more originally "
                "missing entries unimputed."
            )

        return result

    def fit(
        self,
        X: ArrayLike,
        y: Any | None = None,
    ) -> Solver:
        """
        Raise an error because this solver supports only transductive mode.
        """
        del X, y

        raise NotImplementedError(
            f"{self.__class__.__name__} does not support a separate "
            "fit() operation. Use fit_transform() instead."
        )

    def transform(
        self,
        X: ArrayLike,
        y: Any | None = None,
    ) -> FloatArray:
        """
        Raise an error because this solver supports only transductive mode.
        """
        del X, y

        raise NotImplementedError(
            f"{self.__class__.__name__} does not support transform() on "
            "new data. Use fit_transform() instead."
        )

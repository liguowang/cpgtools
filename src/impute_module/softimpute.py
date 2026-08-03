# Portions of this file are derived from "fancyimpute", licensed under
# "Apache License, Version 2.0 "
#
# Modifications Copyright (c) 2024-2026 Liguo Wang
# Modified for use in CpGtools.
#
# See LICENSE.txt and the relevant third-party license notices for details.

from __future__ import annotations

import logging
from numbers import Integral, Real
import numpy as np
from numpy.typing import NDArray
from sklearn.utils import check_array
from sklearn.utils.extmath import randomized_svd
from .common import masked_mae
from .solver import Solver

LOGGER = logging.getLogger(__name__)

FloatArray = NDArray[np.float64]


class SoftImpute(Solver):
    """
    Impute missing matrix entries using nuclear-norm-regularized SVD.

    This implements the SoftImpute procedure described in:

    Mazumder, Hastie, and Tibshirani,
    "Spectral Regularization Algorithms for Learning Large
    Incomplete Matrices."

    At every iteration, observed entries remain fixed while missing entries
    are replaced by values from a soft-thresholded low-rank reconstruction.
    """

    def __init__(
        self,
        shrinkage_value: float | None = None,
        convergence_threshold: float = 1e-3,
        max_iters: int = 100,
        max_rank: int | None = None,
        n_power_iterations: int = 1,
        init_fill_method: str = "zero",
        min_value: float | None = None,
        max_value: float | None = None,
        normalizer: object | None = None,
        random_state: int | None = 1234,
        verbose: bool = True,
    ) -> None:
        """
        Parameters
        ----------
        shrinkage_value
            Nonnegative value subtracted from each singular value. Singular
            values that become nonpositive are removed.

            When ``None``, the value is estimated as one-fiftieth of the
            largest singular value of the initialized matrix.

        convergence_threshold
            Stop when the relative Frobenius-norm change among the originally
            missing entries falls below this value.

        max_iters
            Maximum number of SoftImpute iterations.

        max_rank
            Maximum number of singular components calculated per iteration.
            When supplied, randomized truncated SVD is used. When ``None``,
            a full SVD is calculated.

        n_power_iterations
            Number of power iterations used by randomized SVD.

        init_fill_method
            Initialization method passed to ``Solver``.

        min_value
            Optional lower bound applied only to imputed values. Leave as
            ``None`` for unbounded data such as methylation M values.

        max_value
            Optional upper bound applied only to imputed values. Leave as
            ``None`` for unbounded data such as methylation M values.

        normalizer
            Optional object providing compatible ``fit`` and ``transform``
            methods.

        random_state
            Random seed used by randomized SVD. Use ``None`` for
            nondeterministic decompositions.

        verbose
            Emit per-iteration progress messages through the logging system.
        """
        self._validate_parameters(
            shrinkage_value=shrinkage_value,
            convergence_threshold=convergence_threshold,
            max_iters=max_iters,
            max_rank=max_rank,
            n_power_iterations=n_power_iterations,
            min_value=min_value,
            max_value=max_value,
            random_state=random_state,
            verbose=verbose,
        )

        super().__init__(
            fill_method=init_fill_method,
            min_value=min_value,
            max_value=max_value,
            normalizer=normalizer,
        )

        self.shrinkage_value = shrinkage_value
        self.convergence_threshold = float(convergence_threshold)
        self.max_iters = int(max_iters)
        self.max_rank = None if max_rank is None else int(max_rank)
        self.n_power_iterations = int(n_power_iterations)
        self.random_state = random_state
        self.verbose = verbose

        # Diagnostics populated after solve().
        self.n_iter_: int = 0
        self.rank_: int = 0
        self.shrinkage_value_: float | None = None
        self.converged_: bool = False

    @staticmethod
    def _validate_parameters(
        *,
        shrinkage_value: float | None,
        convergence_threshold: float,
        max_iters: int,
        max_rank: int | None,
        n_power_iterations: int,
        min_value: float | None,
        max_value: float | None,
        random_state: int | None,
        verbose: bool,
    ) -> None:
        """Validate constructor parameters."""
        if shrinkage_value is not None:
            if (
                not isinstance(shrinkage_value, Real)
                or isinstance(shrinkage_value, bool)
                or not np.isfinite(shrinkage_value)
                or shrinkage_value < 0
            ):
                raise ValueError(
                    "shrinkage_value must be a finite, nonnegative number "
                    "or None."
                )

        if (
            not isinstance(convergence_threshold, Real)
            or isinstance(convergence_threshold, bool)
            or not np.isfinite(convergence_threshold)
            or convergence_threshold <= 0
        ):
            raise ValueError(
                "convergence_threshold must be a finite number greater "
                "than zero."
            )

        if (
            not isinstance(max_iters, Integral)
            or isinstance(max_iters, bool)
            or max_iters < 1
        ):
            raise ValueError("max_iters must be a positive integer.")

        if max_rank is not None and (
            not isinstance(max_rank, Integral)
            or isinstance(max_rank, bool)
            or max_rank < 1
        ):
            raise ValueError("max_rank must be a positive integer or None.")

        if (
            not isinstance(n_power_iterations, Integral)
            or isinstance(n_power_iterations, bool)
            or n_power_iterations < 0
        ):
            raise ValueError(
                "n_power_iterations must be a nonnegative integer."
            )

        if min_value is not None and (
            not isinstance(min_value, Real)
            or isinstance(min_value, bool)
            or not np.isfinite(min_value)
        ):
            raise ValueError("min_value must be a finite number or None.")

        if max_value is not None and (
            not isinstance(max_value, Real)
            or isinstance(max_value, bool)
            or not np.isfinite(max_value)
        ):
            raise ValueError("max_value must be a finite number or None.")

        if (
            min_value is not None
            and max_value is not None
            and min_value >= max_value
        ):
            raise ValueError("min_value must be smaller than max_value.")

        if random_state is not None and (
            not isinstance(random_state, Integral)
            or isinstance(random_state, bool)
        ):
            raise TypeError("random_state must be an integer or None.")

        if not isinstance(verbose, bool):
            raise TypeError("verbose must be a boolean.")

    def _relative_change(
        self,
        old: FloatArray,
        new: FloatArray,
        missing_mask: NDArray[np.bool_],
    ) -> float:
        """Calculate normalized change among originally missing entries."""
        old_missing = old[missing_mask]
        new_missing = new[missing_mask]

        if old_missing.size == 0:
            return 0.0

        difference_norm = float(np.linalg.norm(new_missing - old_missing))
        old_norm = float(np.linalg.norm(old_missing))

        denominator = max(old_norm, np.finfo(np.float64).eps)
        return difference_norm / denominator

    def _converged(
        self,
        old: FloatArray,
        new: FloatArray,
        missing_mask: NDArray[np.bool_],
    ) -> tuple[bool, float]:
        """Return convergence state and relative change."""
        relative_change = self._relative_change(old, new, missing_mask)

        return (
            relative_change < self.convergence_threshold,
            relative_change,
        )

    def _clip_imputed_values(
        self,
        values: FloatArray,
    ) -> FloatArray:
        """Apply optional bounds to imputed values only."""
        if self.min_value is None and self.max_value is None:
            return values

        lower = -np.inf if self.min_value is None else self.min_value
        upper = np.inf if self.max_value is None else self.max_value
        return np.clip(values, lower, upper)

    def _svd_step(
        self,
        matrix: FloatArray,
        shrinkage_value: float,
        max_rank: int | None = None,
    ) -> tuple[FloatArray, int]:
        """
        Return a soft-thresholded low-rank reconstruction and its rank.
        """
        min_dimension = min(matrix.shape)

        if max_rank is not None:
            n_components = min(max_rank, min_dimension)

            left, singular_values, right = randomized_svd(
                matrix,
                n_components=n_components,
                n_iter=self.n_power_iterations,
                random_state=self.random_state,
            )
        else:
            left, singular_values, right = np.linalg.svd(
                matrix,
                full_matrices=False,
            )

        thresholded = np.maximum(
            singular_values - shrinkage_value,
            0.0,
        )

        rank = int(np.count_nonzero(thresholded > 0.0))

        if rank == 0:
            return np.zeros_like(matrix), 0

        # Equivalent to U @ diag(s) @ Vt, without constructing diag(s).
        reconstruction = (
            left[:, :rank] * thresholded[:rank]
        ) @ right[:rank, :]

        return reconstruction, rank

    def _max_singular_value(self, matrix: FloatArray) -> float:
        """Estimate the largest singular value of a matrix."""
        _, singular_values, _ = randomized_svd(
            matrix,
            n_components=1,
            n_iter=max(self.n_power_iterations, 5),
            random_state=self.random_state,
        )

        return float(singular_values[0])

    @staticmethod
    def _validate_missing_mask(
        missing_mask: NDArray[np.bool_],
        shape: tuple[int, int],
    ) -> NDArray[np.bool_]:
        """Validate and normalize the missing-value mask."""
        mask = np.asarray(missing_mask)

        if mask.dtype != np.bool_:
            raise TypeError("missing_mask must contain boolean values.")

        if mask.ndim != 2:
            raise ValueError("missing_mask must be two-dimensional.")

        if mask.shape != shape:
            raise ValueError(
                "missing_mask must have the same shape as X: "
                f"expected {shape}, received {mask.shape}."
            )

        return mask

    def solve(
        self,
        X: NDArray[np.floating],
        missing_mask: NDArray[np.bool_],
    ) -> FloatArray:
        """
        Run SoftImpute on an initialized matrix.

        Parameters
        ----------
        X
            Two-dimensional initialized numeric matrix. Missing positions
            should already contain values supplied by the parent ``Solver``.

        missing_mask
            Boolean mask identifying entries that were originally missing.

        Returns
        -------
        numpy.ndarray
            Imputed matrix. Originally observed entries are preserved.
        """
        values = check_array(
            X,
            dtype=np.float64,
            copy=True,
            ensure_2d=True,
            ensure_all_finite=True,
        )

        mask = self._validate_missing_mask(
            missing_mask,
            values.shape,
        )

        if not mask.any():
            self.n_iter_ = 0
            self.rank_ = int(np.linalg.matrix_rank(values))
            self.shrinkage_value_ = self.shrinkage_value
            self.converged_ = True
            return values

        min_dimension = min(values.shape)

        if self.max_rank is not None and self.max_rank > min_dimension:
            raise ValueError(
                "max_rank cannot exceed the smaller matrix dimension "
                f"({min_dimension}); received {self.max_rank}."
            )

        # Keep an immutable copy for observed-value diagnostics.
        initialized = values.copy()

        # This copy is updated only at the originally missing positions.
        filled = values.copy()
        observed_mask = ~mask

        maximum_singular_value = self._max_singular_value(filled)

        if self.shrinkage_value is None:
            shrinkage_value = maximum_singular_value / 50.0
        else:
            shrinkage_value = float(self.shrinkage_value)

        self.shrinkage_value_ = shrinkage_value
        self.converged_ = False
        self.n_iter_ = 0
        self.rank_ = 0

        if self.verbose:
            LOGGER.info(
                "SoftImpute maximum initial singular value: %.6g",
                maximum_singular_value,
            )
            LOGGER.info(
                "SoftImpute shrinkage value: %.6g",
                shrinkage_value,
            )

        for iteration in range(1, self.max_iters + 1):
            reconstruction, rank = self._svd_step(
                filled,
                shrinkage_value,
                max_rank=self.max_rank,
            )

            # Apply optional bounds only to entries that will be imputed.
            # Observed values and their low-rank reconstructions are not
            # clipped. With both limits set to None, M-value predictions
            # remain unbounded.
            updated = filled.copy()
            updated[mask] = self._clip_imputed_values(
                reconstruction[mask]
            )

            converged, relative_change = self._converged(
                old=filled,
                new=updated,
                missing_mask=mask,
            )

            if self.verbose:
                observed_mae = masked_mae(
                    X_true=initialized,
                    X_pred=reconstruction,
                    mask=observed_mask,
                )

                LOGGER.info(
                    "SoftImpute iteration %d: observed MAE=%.6g, "
                    "relative change=%.6g, rank=%d",
                    iteration,
                    observed_mae,
                    relative_change,
                    rank,
                )

            # Only missing entries are updated. Observed values remain fixed.
            filled = updated

            self.n_iter_ = iteration
            self.rank_ = rank

            if converged:
                self.converged_ = True
                break

        if self.verbose:
            if self.converged_:
                LOGGER.info(
                    "SoftImpute converged after %d iteration(s), "
                    "lambda=%.6g, rank=%d.",
                    self.n_iter_,
                    shrinkage_value,
                    self.rank_,
                )
            else:
                LOGGER.warning(
                    "SoftImpute reached %d iteration(s) without meeting "
                    "the convergence threshold; returning the final "
                    "estimates. Lambda=%.6g, rank=%d.",
                    self.max_iters,
                    shrinkage_value,
                    self.rank_,
                )

        if not np.isfinite(filled).all():
            raise RuntimeError(
                "SoftImpute produced non-finite values."
            )

        return filled

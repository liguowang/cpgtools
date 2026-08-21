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
from collections import defaultdict
from numbers import Integral
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


LOGGER = logging.getLogger(__name__)

GKNNMethod = Literal["AA", "WA", "TA"]
CpGAnnotation = tuple[str, int, int, frozenset[str]]
GenomicIndex = dict[str, tuple[np.ndarray, np.ndarray]]


def _load_cpg_annotations(
    gfile: str | Path,
) -> tuple[dict[str, CpGAnnotation], GenomicIndex]:
    """Read CpG annotations and build chromosome-specific coordinate indexes."""
    path = Path(gfile)
    annotations: dict[str, CpGAnnotation] = {}
    chromosome_records: dict[str, list[tuple[int, str]]] = defaultdict(list)

    try:
        with path.open(encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()

                if not line or line.startswith("#"):
                    continue

                fields = line.split()
                if len(fields) < 5:
                    LOGGER.warning(
                        "Skipping line %d of %s: fewer than five columns.",
                        line_number,
                        path,
                    )
                    continue

                chrom, start_text, end_text, cpg_id, cre_text = fields[:5]

                try:
                    start = int(start_text)
                    end = int(end_text)
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid genomic coordinate on line {line_number} "
                        f"of {path}."
                    ) from exc

                if start < 0 or end < start:
                    raise ValueError(
                        f"Invalid genomic interval on line {line_number} "
                        f"of {path}: {chrom}:{start}-{end}."
                    )

                if cpg_id in annotations:
                    raise ValueError(
                        f"CpG {cpg_id!r} occurs more than once in {path}."
                    )

                if cre_text.upper() in {"N/A", "NA", ".", "NONE"}:
                    cres = frozenset()
                else:
                    cres = frozenset(
                        cre.strip()
                        for cre in cre_text.split(",")
                        if cre.strip()
                    )

                annotations[cpg_id] = (chrom, start, end, cres)
                chromosome_records[chrom].append((start, cpg_id))
    except OSError as exc:
        raise ValueError(f"Unable to read genomic annotation file {path}: {exc}") from exc

    if not annotations:
        raise ValueError(f"No valid CpG annotations were read from {path}.")

    genomic_index: GenomicIndex = {}

    for chrom, records in chromosome_records.items():
        records.sort(key=lambda item: (item[0], item[1]))
        starts = np.asarray([start for start, _ in records], dtype=np.int64)
        cpg_ids = np.asarray([cpg_id for _, cpg_id in records], dtype=object)
        genomic_index[chrom] = (starts, cpg_ids)

    return annotations, genomic_index


def _read_new_cpgs(cpgfile: str | Path) -> list[str]:
    """Read additional CpG IDs, ignoring comments, blanks, and duplicates."""
    path = Path(cpgfile)
    cpg_ids: list[str] = []
    seen: set[str] = set()

    try:
        with path.open(encoding="utf-8") as handle:
            for raw_line in handle:
                cpg_id = raw_line.strip()
                if not cpg_id or cpg_id.startswith("#") or cpg_id in seen:
                    continue
                cpg_ids.append(cpg_id)
                seen.add(cpg_id)
    except OSError as exc:
        raise ValueError(f"Unable to read CpG file {path}: {exc}") from exc

    return cpg_ids


def _find_genomic_neighbors(
    target_cpg: str,
    annotations: dict[str, CpGAnnotation],
    genomic_index: GenomicIndex,
    available_cpgs: set[str],
    *,
    up_dist: int,
    down_dist: int,
    up_ncpg: int,
    down_ncpg: int,
    same_cre: bool,
) -> list[tuple[str, int]]:
    """Find the closest eligible upstream and downstream CpGs."""
    target_annotation = annotations.get(target_cpg)
    if target_annotation is None:
        return []

    target_chrom, target_start, _, target_cres = target_annotation

    if same_cre and not target_cres:
        return []

    chromosome_index = genomic_index.get(target_chrom)
    if chromosome_index is None:
        return []

    starts, cpg_ids = chromosome_index

    left = int(np.searchsorted(starts, target_start - up_dist, side="left"))
    right = int(np.searchsorted(starts, target_start + down_dist, side="right"))

    upstream: list[tuple[int, str]] = []
    downstream: list[tuple[int, str]] = []

    for candidate_start, candidate_cpg in zip(
        starts[left:right],
        cpg_ids[left:right],
    ):
        candidate_cpg = str(candidate_cpg)

        if candidate_cpg == target_cpg or candidate_cpg not in available_cpgs:
            continue

        candidate_cres = annotations[candidate_cpg][3]
        if same_cre and not target_cres.intersection(candidate_cres):
            continue

        difference = int(candidate_start) - target_start

        if difference < 0:
            upstream.append((-difference, candidate_cpg))
        elif difference > 0:
            downstream.append((difference, candidate_cpg))
        else:
            # Distinct CpGs sharing a start coordinate are treated as one
            # base pair away to avoid zero-distance weights.
            downstream.append((1, candidate_cpg))

    upstream.sort(key=lambda item: (item[0], item[1]))
    downstream.sort(key=lambda item: (item[0], item[1]))

    selected = upstream[:up_ncpg] + downstream[:down_ncpg]
    return [(cpg_id, distance) for distance, cpg_id in selected]


def _weighted_average(values: np.ndarray, distances: np.ndarray) -> float:
    """Calculate an inverse-distance weighted average."""
    if values.size == 0:
        return float("nan")

    safe_distances = np.maximum(distances.astype(float), 1.0)
    weights = 1.0 / safe_distances
    return float(np.average(values, weights=weights))


def _trimmed_average(values: np.ndarray) -> float:
    """Remove values farther than two population SDs from the mean."""
    if values.size == 0:
        return float("nan")

    if values.size < 3:
        return float(np.mean(values))

    mean = float(np.mean(values))
    standard_deviation = float(np.std(values, ddof=0))

    if standard_deviation == 0.0:
        return mean

    retained = values[np.abs(values - mean) <= 2.0 * standard_deviation]
    if retained.size == 0:
        return mean

    return float(np.mean(retained))


def fill_GNN(
    df: pd.DataFrame,
    gfile: str | Path,
    cpgfile: str | Path | None = None,
    *,
    up_dist: int = 100,
    down_dist: int = 100,
    up_ncpg: int = 2,
    down_ncpg: int = 2,
    same_CRE: bool = False,
    method: GKNNMethod = "TA",) -> pd.DataFrame:
    """
    Impute missing methylation values using neighboring genomic CpGs.

    CpGs are searched independently on each chromosome using sorted genomic
    coordinates and binary search. For each target CpG, at most ``up_ncpg``
    upstream and ``down_ncpg`` downstream CpGs within the requested distance
    limits are used. When ``same_CRE`` is true, neighbors must share at least
    one candidate regulatory element with the target CpG.

    Parameters
    ----------
    df
        DNA methylation matrix with CpGs in rows and samples in columns.

    gfile
        Genomic annotation file with at least five whitespace-delimited
        columns: ``chrom``, ``start``, ``end``, ``cpg_id``, and ``CRE``.
        Multiple CREs must be comma-separated. Use ``N/A`` for no CRE.

    cpgfile
        Optional file containing additional CpG IDs, one per line. These CpGs
        are appended as missing rows and imputed when possible.

    up_dist, down_dist
        Maximum upstream and downstream search distances in base pairs.

    up_ncpg, down_ncpg
        Maximum numbers of upstream and downstream neighbors.

    same_CRE
        If true, neighbors must share at least one CRE with the target CpG.

    method
        ``"AA"`` for arithmetic average, ``"WA"`` for inverse-distance
        weighted average, or ``"TA"`` for a two-standard-deviation trimmed
        average.

    Returns
    -------
    pandas.DataFrame
        A copy of the matrix with imputable values filled. Rows that cannot
        be imputed remain missing.

    Notes
    -----
    Only values present before imputation are used as donors. Newly imputed
    values are not reused for other targets, so results do not depend on row
    processing order.
    """
    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    if df.empty:
        raise ValueError("df must not be empty.")

    if df.index.has_duplicates:
        raise ValueError("df contains duplicated CpG identifiers.")

    if df.columns.has_duplicates:
        raise ValueError("df contains duplicated sample names.")

    if method not in {"AA", "WA", "TA"}:
        raise ValueError("method must be one of 'AA', 'WA', or 'TA'.")

    integer_parameters = {
        "up_dist": up_dist,
        "down_dist": down_dist,
        "up_ncpg": up_ncpg,
        "down_ncpg": down_ncpg,
    }
    for name, value in integer_parameters.items():
        if not isinstance(value, Integral) or isinstance(value, bool):
            raise TypeError(f"{name} must be an integer.")
        if value < 0:
            raise ValueError(f"{name} must be nonnegative.")

    if not isinstance(same_CRE, bool):
        raise TypeError("same_CRE must be a Boolean value.")

    try:
        working_df = df.astype(float).copy()
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "df must contain only numeric values and missing values."
        ) from exc

    if np.isinf(working_df.to_numpy()).any():
        raise ValueError("df must not contain infinite values.")

    if cpgfile is not None:
        requested_cpgs = _read_new_cpgs(cpgfile)
        existing = set(working_df.index)
        new_cpgs = [cpg_id for cpg_id in requested_cpgs if cpg_id not in existing]

        if new_cpgs:
            LOGGER.info("Adding %d CpGs from %s.", len(new_cpgs), cpgfile)
            working_df = working_df.reindex([*working_df.index, *new_cpgs])

    annotations, genomic_index = _load_cpg_annotations(gfile)

    # Freeze donor values before any imputation to prevent chained imputation.
    donor_df = working_df.copy()
    result = working_df.copy()
    available_cpgs = set(donor_df.index)

    target_cpgs = result.index[result.isna().any(axis=1)].tolist()
    LOGGER.info(
        "Searching for genomic neighbors for %d CpGs with missing values.",
        len(target_cpgs),
    )

    neighbor_cache: dict[str, list[tuple[str, int]]] = {}
    missing_annotation_count = 0

    for target_cpg in target_cpgs:
        if target_cpg not in annotations:
            missing_annotation_count += 1

        neighbor_cache[target_cpg] = _find_genomic_neighbors(
            target_cpg=target_cpg,
            annotations=annotations,
            genomic_index=genomic_index,
            available_cpgs=available_cpgs,
            up_dist=up_dist,
            down_dist=down_dist,
            up_ncpg=up_ncpg,
            down_ncpg=down_ncpg,
            same_cre=same_CRE,
        )

    if missing_annotation_count:
        LOGGER.warning(
            "%d target CpGs are absent from the annotation file.",
            missing_annotation_count,
        )

    imputed_count = 0

    for target_cpg in target_cpgs:
        neighbors = neighbor_cache[target_cpg]
        if not neighbors:
            continue

        missing_samples = result.columns[result.loc[target_cpg].isna()]

        for sample in missing_samples:
            values: list[float] = []
            distances: list[float] = []

            for neighbor_cpg, distance in neighbors:
                value = donor_df.at[neighbor_cpg, sample]
                if pd.isna(value):
                    continue

                values.append(float(value))
                distances.append(float(distance))

            if not values:
                continue

            value_array = np.asarray(values, dtype=float)
            distance_array = np.asarray(distances, dtype=float)

            if method == "AA":
                imputed_value = float(np.mean(value_array))
            elif method == "WA":
                imputed_value = _weighted_average(value_array, distance_array)
            else:
                imputed_value = _trimmed_average(value_array)

            result.at[target_cpg, sample] = imputed_value
            imputed_count += 1

    LOGGER.info(
        "GNN imputed %d values; %d missing values remain in %d "
        "entirely missing rows.",
        imputed_count,
        int(result.isna().to_numpy().sum()),
        int(result.isna().all(axis=1).sum()),
    )

    return result


# Optional lowercase alias following PEP 8 naming conventions.
fill_gnn = fill_GNN


__all__ = ["fill_gnn", "fill_GNN"]

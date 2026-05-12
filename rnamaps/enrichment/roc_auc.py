"""ROC/AUC analysis of CLIP signal as a predictor of exon regulation.

For each splice-site region and each non-control category, we treat
the per-exon CLIP signal as a continuous score and ask how well it
discriminates the category from the control set. Two complementary
outputs are produced:

1. **Per-position AUC line** (``plot_kind="roc_auc"``, plot_df). At
   every base ``p`` along the splice-site window we use the per-exon
   signal at that single base as the prediction score and the binary
   {category, control} membership as the label, then compute the
   rank-based AUC (equivalent to the Mann–Whitney U statistic divided
   by the product of group sizes). Plotted as a line vs position with
   the no-information baseline at 0.5 (or, equivalently, the signed
   ``2*(AUC - 0.5)`` transform centred on 0). This answers *where*
   around the splice site the score discriminates regulated exons.

2. **Per-region ROC curve** (``extras['roc_curves_df']``). The per-exon
   signal is aggregated over the window (``mean`` and ``max`` by
   default) into a single score per exon, then a standard ROC curve
   (TPR vs FPR over thresholds) and its AUC are computed across exons
   in ``{category, control}``. Plotted as one panel per
   (region, category, aggregator) with the AUC printed in the legend.

A label-permutation p-value can optionally be produced for both
outputs via ``n_perm > 0`` -- shuffles the ``{category, control}``
labels ``n_perm`` times and reports the proportion of permutations
with an AUC at least as far from 0.5 as the observed one (two-sided).

The same per-exon coverage matrix used by every other enrichment
method is consumed here; if the CLIP BED carries meaningful scores,
they need to be plumbed into the matrix at the coverage step via
``--xl_score`` -- this module does not re-read the BED.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats as _sstats

from rnamaps.coverage import aggregate_legacy_columns
from rnamaps.enrichment._base import EnrichmentResult
from rnamaps.permutation import _build_coverage_matrix


def _rank_auc(scores: np.ndarray, labels: np.ndarray) -> float:
    """AUC of ``scores`` for binary ``labels`` (1=positive, 0=negative).

    Computed via the Mann–Whitney U identity ``AUC = U / (n_pos * n_neg)``
    using ``scipy.stats.rankdata`` for tie handling (mid-ranks):

        AUC = (R_pos - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

    where ``R_pos`` is the sum of ranks of positive-class scores in the
    pooled ranking. Returns 0.5 when either class is empty so callers
    don't have to special-case degenerate slices.
    """
    n = labels.size
    n_pos = int(labels.sum())
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    ranks = _sstats.rankdata(scores, method="average")
    r_pos = ranks[labels.astype(bool)].sum()
    return float(
        (r_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    )


def _vectorised_auc_per_column(matrix: np.ndarray, is_pos: np.ndarray) -> np.ndarray:
    """Per-column AUC for a ``(n_exons, n_positions)`` score matrix.

    Vectorised over positions: ranks each column independently and
    applies the Mann–Whitney formula across all columns at once.

    Parameters
    ----------
    matrix : ndarray, shape (n_exons, n_positions)
        Per-exon score (e.g. the coverage matrix). ``float`` recommended
        for stable rank tie behaviour with non-binarised scores.
    is_pos : ndarray of bool, shape (n_exons,)
        Positive-class membership (e.g. ``exon_names == 'enhanced'``).

    Returns
    -------
    ndarray, shape (n_positions,)
        AUC per position; positions where either class is empty get 0.5.
    """
    n_pos = int(is_pos.sum())
    n_neg = int((~is_pos).sum())
    if n_pos == 0 or n_neg == 0 or matrix.size == 0:
        return np.full(matrix.shape[1], 0.5, dtype=np.float64)
    ranks = _sstats.rankdata(matrix, method="average", axis=0)
    r_pos = ranks[is_pos].sum(axis=0)
    auc = (r_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return np.asarray(auc, dtype=np.float64)


def _roc_curve(scores: np.ndarray, labels: np.ndarray):
    """Threshold-sweep ROC curve points without a sklearn dependency.

    Returns ``(fpr, tpr)`` sorted by decreasing threshold, with the
    canonical anchor points (0,0) prepended and (1,1) appended so the
    curve renders to the full unit square.
    """
    n_pos = int(labels.sum())
    n_neg = int(labels.size - n_pos)
    if n_pos == 0 or n_neg == 0 or scores.size == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])
    order = np.argsort(-scores, kind="mergesort")
    s = scores[order]
    y = labels[order].astype(np.int64, copy=False)
    # Cumulative TP / FP as the threshold drops past each score.
    tps = np.cumsum(y)
    fps = np.cumsum(1 - y)
    # Collapse ties: only keep the last point at each unique threshold,
    # otherwise vertical/horizontal segments are drawn as staircases.
    distinct_mask = np.r_[np.diff(s) != 0, True]
    tps = tps[distinct_mask]
    fps = fps[distinct_mask]
    tpr = np.r_[0.0, tps / n_pos]
    fpr = np.r_[0.0, fps / n_neg]
    # Append (1,1) if not already there (mergesort-tied tail).
    if tpr[-1] < 1.0 or fpr[-1] < 1.0:
        tpr = np.r_[tpr, 1.0]
        fpr = np.r_[fpr, 1.0]
    return fpr, tpr


def _aggregate_window(matrix: np.ndarray, how: str) -> np.ndarray:
    """Collapse a ``(n_exons, n_positions)`` matrix to a per-exon score."""
    if how == "mean":
        return matrix.mean(axis=1)
    if how == "max":
        return matrix.max(axis=1)
    if how == "sum":
        return matrix.sum(axis=1)
    raise ValueError(f"Unknown ROC aggregator: {how!r}")


def _resolve_aggregators(roc_aggregator: str):
    if roc_aggregator == "both":
        return ("mean", "max")
    if roc_aggregator in ("mean", "max", "sum"):
        return (roc_aggregator,)
    raise ValueError(
        f"roc_aggregator must be one of {{'mean','max','both'}}, got "
        f"{roc_aggregator!r}"
    )


def _permutation_auc_pvalue(
    matrix: np.ndarray,
    is_cat: np.ndarray,
    observed: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Two-sided permutation p-value for per-position AUC.

    Shuffles the {category, control} labels ``n_perm`` times and
    computes ``2 * min(P(null >= |obs - 0.5|), 1)`` per position via the
    pivot ``|null_auc - 0.5|``. Vectorised across positions.
    """
    if n_perm <= 0:
        return np.full(matrix.shape[1], np.nan, dtype=np.float64)
    n_total = matrix.shape[0]
    n_c = int(is_cat.sum())
    obs_dev = np.abs(observed - 0.5)
    null_ge = np.zeros(matrix.shape[1], dtype=np.int64)
    for _ in range(n_perm):
        idx = rng.choice(n_total, size=n_c, replace=False)
        perm_mask = np.zeros(n_total, dtype=bool)
        perm_mask[idx] = True
        null_auc = _vectorised_auc_per_column(matrix, perm_mask)
        null_ge += (np.abs(null_auc - 0.5) >= obs_dev).astype(np.int64)
    return (null_ge + 1) / (n_perm + 1)


def compute(
    region,
    exon_categories: pd.Series,
    label: str,
    *,
    rng: np.random.Generator,
    roc_aggregator: str = "both",
    n_perm: int = 0,
    smoothing: int = 1,
    control_label: str = "control",
    binarise: bool = False,
) -> EnrichmentResult:
    """Compute per-position AUC and per-region ROC for one splice-site region.

    Parameters
    ----------
    region : rnamaps.coverage.RegionCoverage or pd.DataFrame
        Per-exon coverage / signal matrix for one splice-site region.
    exon_categories : Series
        Counts per category. Must include ``control_label``.
    label : str
        Splice-site region label (e.g. ``"middle_3ss"``).
    rng : np.random.Generator
    roc_aggregator : {"mean", "max", "both"}
        Per-exon aggregation for the per-region ROC curve.
    n_perm : int
        Number of label permutations for the AUC null. ``0`` (default)
        skips the null and the ``pvalue`` columns are NaN.
    smoothing : int
        Centred Gaussian rolling-mean window applied to the per-position
        AUC line (matches the smoothing used elsewhere). Set ``<=1`` to
        disable.
    control_label : str
    binarise : bool
        ROC/AUC is rank-based and invariant to any monotone rescaling,
        so binarisation usually does not help here. Defaults to
        ``False`` (use raw signal) regardless of the global
        ``--binarise`` setting; pass ``True`` for a sanity check of
        "presence-only" classification.

    Returns
    -------
    EnrichmentResult
        ``plot_df`` columns: ``name, position, label, auc,
        auc_signed, pvalue, coverage, ...``. ``auc_signed`` is
        ``2 * (auc - 0.5)`` so the no-info baseline is 0 and the sign
        matches the rest of the RNA maps.
        ``extras['roc_curves_df']`` carries one row per
        (region, category, aggregator, threshold-step) for plotting.
        ``extras['region_auc_df']`` carries one row per
        (region, category, aggregator) with the headline AUC.
    """
    if control_label not in exon_categories.index:
        raise ValueError(
            f"roc_auc requires a '{control_label}' category."
        )

    categories = [c for c in exon_categories.index if c != control_label]
    n_ctrl = int(exon_categories.loc[control_label])
    aggregators = _resolve_aggregators(roc_aggregator)

    agg = aggregate_legacy_columns(region, exon_categories, control_label)

    plot_rows = []
    roc_rows = []
    region_auc_rows = []

    for cat in categories:
        matrix, is_cat, positions = _build_coverage_matrix(
            region, cat, control_label, binarise=binarise,
        )
        if matrix is None or is_cat.sum() == 0 or n_ctrl == 0:
            logging.warning(
                f"[roc_auc] Skipping {cat} ({label}): empty matrix."
            )
            continue

        # Per-position AUC (over the {category, control} subset).
        auc_pos = _vectorised_auc_per_column(matrix, is_cat)
        auc_signed = 2.0 * (auc_pos - 0.5)
        if smoothing and smoothing > 1:
            auc_signed_smoothed = pd.Series(auc_signed).rolling(
                smoothing, center=True, win_type='gaussian',
            ).mean(std=2).to_numpy()
        else:
            auc_signed_smoothed = auc_signed.copy()
        if n_perm > 0:
            pvals = _permutation_auc_pvalue(
                matrix, is_cat, auc_pos, n_perm, rng,
            )
        else:
            pvals = np.full(matrix.shape[1], np.nan, dtype=np.float64)

        plot_rows.append(pd.DataFrame({
            'name': cat,
            'position': positions,
            'label': label,
            'auc': auc_pos,
            'auc_signed': auc_signed,
            'auc_signed_smoothed': auc_signed_smoothed,
            'pvalue': pvals,
        }))

        # Per-region ROC.
        labels_bin = is_cat.astype(np.int64, copy=False)
        for how in aggregators:
            scores = _aggregate_window(matrix, how)
            auc_region = _rank_auc(scores, labels_bin)
            fpr, tpr = _roc_curve(scores, labels_bin)
            # Optional permutation null on the region AUC (cheap; uses
            # the same shuffles philosophy as the per-position path).
            if n_perm > 0:
                obs_dev = abs(auc_region - 0.5)
                ge = 0
                n_total = matrix.shape[0]
                n_c = int(is_cat.sum())
                for _ in range(n_perm):
                    idx = rng.choice(n_total, size=n_c, replace=False)
                    mask = np.zeros(n_total, dtype=bool)
                    mask[idx] = True
                    null_auc = _rank_auc(
                        scores, mask.astype(np.int64, copy=False)
                    )
                    if abs(null_auc - 0.5) >= obs_dev:
                        ge += 1
                p_region = (ge + 1) / (n_perm + 1)
            else:
                p_region = np.nan

            roc_rows.append(pd.DataFrame({
                'name': cat,
                'label': label,
                'aggregator': how,
                'fpr': fpr,
                'tpr': tpr,
            }))
            region_auc_rows.append({
                'name': cat,
                'label': label,
                'aggregator': how,
                'auc': float(auc_region),
                'n_cat': int(is_cat.sum()),
                'n_ctrl': int((~is_cat).sum()),
                'pvalue': float(p_region) if p_region == p_region
                          else np.nan,
            })

    # Control rows on the per-position line plot (zero baseline for
    # legend continuity).
    ctrl_positions = agg.loc[agg['name'] == control_label, 'position'].to_numpy()
    if len(ctrl_positions):
        plot_rows.append(pd.DataFrame({
            'name': control_label,
            'position': ctrl_positions,
            'label': label,
            'auc': 0.5,
            'auc_signed': 0.0,
            'auc_signed_smoothed': 0.0,
            'pvalue': np.nan,
        }))

    if not plot_rows:
        plot_df = pd.DataFrame(columns=[
            'name', 'position', 'label', 'auc', 'auc_signed',
            'auc_signed_smoothed', 'pvalue',
        ])
    else:
        plot_df = pd.concat(plot_rows, ignore_index=True)

    plot_df = plot_df.merge(
        agg[['name', 'position', 'coverage', 'number_exons', 'norm_coverage',
             'control_coverage', 'control_number_exons',
             'control_norm_coverage', 'fold_change']],
        on=['name', 'position'], how='left'
    )

    roc_curves_df = (pd.concat(roc_rows, ignore_index=True)
                     if roc_rows else pd.DataFrame(columns=[
                         'name', 'label', 'aggregator', 'fpr', 'tpr']))
    region_auc_df = (pd.DataFrame(region_auc_rows)
                     if region_auc_rows else pd.DataFrame(columns=[
                         'name', 'label', 'aggregator', 'auc',
                         'n_cat', 'n_ctrl', 'pvalue']))

    return EnrichmentResult(
        plot_df=plot_df,
        clusters_df=pd.DataFrame(),
        plot_kind='roc_auc',
        y_columns=('auc_signed_smoothed',),
        method_name='roc_auc',
        ylabel='signed AUC (2*(AUC - 0.5)) vs control',
        extras={
            'roc_curves_df': roc_curves_df,
            'region_auc_df': region_auc_df,
        },
    )

"""Label-permutation z-score test for position-wise CLIP enrichment.

Replaces the per-position Fisher's exact test against control with a
standardised score derived from randomly relabelling exons between the
focal category (e.g. ``enhanced``) and ``control``. For each position
the observed mean-coverage difference is converted to a z-score against
the permutation null mean and standard deviation, and a two-sided
p-value is reported via a normal-tail approximation. This keeps the
reported -log10(p) continuous and unbounded so it can be compared
across datasets and peak callers.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
from scipy import stats as _sstats

from rnamaps.coverage import aggregate_legacy_columns as _aggregate_legacy_columns


# Tiny floor to avoid log10(0) at extreme z-scores.
_PVAL_FLOOR = 1e-300


def _build_coverage_matrix(region, category: str,
                           control_label: str = "control",
                           binarise: bool = True):
    """Select ``{category, control}`` rows from a :class:`RegionCoverage`.

    Parameters
    ----------
    region : rnamaps.coverage.RegionCoverage or pd.DataFrame
        Either the new compact value object or, for backward compat with
        existing tests, the legacy long-form per-exon DataFrame.
    category : str
        Non-control category to compare against ``control_label``.
    control_label : str
    binarise : bool, default True
        When ``True`` (the historical default), the matrix is thresholded
        to 0/1 so downstream methods operate on "did this exon have any
        signal at this position?" semantics. This is the right choice
        for CLIP crosslinks, where raw counts at a single base are
        dominated by sequencing-depth noise and the meaningful signal is
        whether or not an exon was bound. When ``False``, the raw values
        from the coverage matrix are passed through as ``float64``,
        which is what you want for continuous inputs (e.g. AI prediction
        scores, density tracks, or any matrix where the per-(exon,
        position) magnitude is itself informative). The downstream test
        statistic ``mean(cat) - mean(ctrl)`` is then literally a
        difference of mean signal levels rather than a difference of
        positive-exon fractions.

    Returns
    -------
    matrix : np.ndarray, shape (n_rows, n_positions)
        ``float64`` coverage. Rows in ``{category, control_label}``.
        Binarised iff ``binarise`` is ``True``.
    is_category : np.ndarray (bool)
    positions : np.ndarray
    """
    # Backward-compat: accept the legacy long-form DataFrame too.
    if isinstance(region, pd.DataFrame):
        sub = region[region['name'].isin([category, control_label])]
        if sub.empty:
            return None, None, None
        pivot = sub.pivot_table(
            index=['name', 'exon_id'], columns='position', values='coverage',
            fill_value=0, aggfunc='sum'
        ).sort_index(axis=1)
        matrix = pivot.to_numpy(dtype=np.float64)
        if binarise:
            matrix = (matrix > 0).astype(np.float64, copy=False)
        is_category = np.asarray(
            pivot.index.get_level_values('name') == category
        )
        positions = pivot.columns.to_numpy()
        return matrix, is_category, positions

    mask = (region.exon_names == category) | (region.exon_names == control_label)
    if not mask.any():
        return None, None, None
    sub_mat = region.matrix[mask]
    if sub_mat.size == 0:
        return None, None, None
    if binarise:
        matrix = (sub_mat > 0).astype(np.float64, copy=False)
    else:
        matrix = sub_mat.astype(np.float64, copy=False)
    sub_names = region.exon_names[mask]
    is_category = (sub_names == category)
    return matrix, is_category, region.positions.copy()


def _permutation_null(matrix: np.ndarray, is_category: np.ndarray,
                      n_perm: int, rng: np.random.Generator):
    """Run a vectorised label-permutation test on a coverage matrix.

    Returns
    -------
    t_obs : np.ndarray, shape (n_positions,)
    t_null : np.ndarray, shape (n_perm, n_positions)
        Per-permutation, per-position test statistic under the null.
    """
    n_total, n_pos = matrix.shape
    n_cat = int(is_category.sum())
    n_ctrl = n_total - n_cat

    sum_all = matrix.sum(axis=0)
    sum_cat_obs = matrix[is_category].sum(axis=0)
    t_obs = sum_cat_obs / n_cat - (sum_all - sum_cat_obs) / n_ctrl

    t_null = np.empty((n_perm, n_pos), dtype=np.float64)
    for b in range(n_perm):
        idx = rng.choice(n_total, size=n_cat, replace=False)
        sum_cat_perm = matrix[idx].sum(axis=0)
        t_null[b] = sum_cat_perm / n_cat - (sum_all - sum_cat_perm) / n_ctrl

    return t_obs, t_null


def _zscore_two_sided_p(t_obs: np.ndarray, t_null: np.ndarray):
    """Per-position z-score against the permutation null and two-sided p.

    The z-score is computed as ``(t_obs - mean(t_null)) / sd(t_null)`` and
    converted to a two-sided p-value via a standard-normal tail. The
    permutation null for ``mean(coverage_cat) - mean(coverage_ctrl)`` is
    well-approximated by a normal for realistic exon counts (CLT applies
    position-wise across exons), so this gives a continuous, unbounded
    score that does not saturate at ``1 / (B + 1)``.

    Returns
    -------
    z : np.ndarray, shape (n_positions,)
    pvalues : np.ndarray, shape (n_positions,)
    """
    mu = t_null.mean(axis=0)
    sd = t_null.std(axis=0, ddof=1)
    # Guard positions where the null is degenerate (e.g. zero coverage).
    z = np.zeros_like(t_obs, dtype=np.float64)
    valid = sd > 0
    z[valid] = (t_obs[valid] - mu[valid]) / sd[valid]
    pvalues = 2.0 * _sstats.norm.sf(np.abs(z))
    return z, pvalues


def _smooth(values: np.ndarray, smoothing: int) -> np.ndarray:
    if smoothing is None or smoothing <= 1:
        return values.copy()
    s = pd.Series(values)
    return s.rolling(smoothing, center=True, win_type='gaussian').mean(std=2).to_numpy()


def compute_permutation_pvalues(region,
                                exon_categories: pd.Series,
                                label: str,
                                n_perm: int,
                                smoothing: int,
                                rng: np.random.Generator,
                                control_label: str = "control",
                                binarise: bool = True):
    """Compute permutation p-values for every non-control category in a region.

    Parameters
    ----------
    region : rnamaps.coverage.RegionCoverage or pd.DataFrame
        Per-exon coverage for one splice-site region. Long-form
        DataFrame input is still accepted for backward compatibility
        with existing tests.
    exon_categories : Series
        Counts per category (used to detect available categories).
    label : str
        Splice-site region label (e.g. ``"middle_3ss"``).
    n_perm : int
        Number of permutations.
    smoothing : int
        Gaussian rolling window applied to signed -log10(p).
    rng : np.random.Generator
    control_label : str
    binarise : bool, default True
        Whether to threshold the per-exon coverage matrix to 0/1 before
        running the test (CLIP-style "did this exon have any signal at
        this base?" semantics). Set to ``False`` for continuous inputs
        such as AI prediction scores where the per-position magnitude
        is itself meaningful; the test statistic then becomes a
        difference of mean signal levels rather than a difference of
        positive-exon fractions.

    Returns
    -------
    plot_df : DataFrame
        Columns: name, position, label, coverage, number_exons, norm_coverage,
        control_coverage, control_number_exons, control_norm_coverage,
        fold_change, T_obs, zscore, pvalue, -log10pvalue, -log10pvalue_smoothed.
        Includes a row per (category, position) for non-control categories
        plus the control rows themselves (for legend / line continuity).
    clusters_df : DataFrame
        Always empty; returned for backward compatibility with the previous
        signature.
    """
    if control_label not in exon_categories.index:
        raise ValueError(
            f"Permutation test requires a '{control_label}' category."
        )

    categories = [c for c in exon_categories.index if c != control_label]
    n_ctrl = int(exon_categories.loc[control_label])

    agg = _aggregate_legacy_columns(region, exon_categories, control_label)

    plot_rows = []

    for cat in categories:
        matrix, is_cat, positions = _build_coverage_matrix(
            region, cat, control_label, binarise=binarise,
        )
        if matrix is None or is_cat.sum() == 0 or n_ctrl == 0:
            logging.warning(
                f"[permutation] Skipping {cat} ({label}): empty matrix.")
            continue

        t_obs, t_null = _permutation_null(matrix, is_cat, n_perm, rng)
        z, pvalues = _zscore_two_sided_p(t_obs, t_null)

        signed_log10p = -np.log10(np.maximum(pvalues, _PVAL_FLOOR))
        # Sign by direction of z (positive = enriched over control).
        signed_log10p = np.where(z >= 0, signed_log10p, -signed_log10p)

        smoothed = _smooth(signed_log10p, smoothing)

        cat_df = pd.DataFrame({
            'name': cat,
            'position': positions,
            'label': label,
            'T_obs': t_obs,
            'zscore': z,
            'pvalue': pvalues,
            '-log10pvalue': signed_log10p,
            '-log10pvalue_smoothed': smoothed,
        })
        plot_rows.append(cat_df)

    # Add control rows for legend continuity (zeros on the y-axis).
    ctrl_positions = agg.loc[agg['name'] == control_label, 'position'].to_numpy()
    if len(ctrl_positions):
        plot_rows.append(pd.DataFrame({
            'name': control_label,
            'position': ctrl_positions,
            'label': label,
            'T_obs': 0.0,
            'zscore': 0.0,
            'pvalue': 1.0,
            '-log10pvalue': 0.0,
            '-log10pvalue_smoothed': 0.0,
        }))

    if not plot_rows:
        plot_df = pd.DataFrame(columns=[
            'name', 'position', 'label', 'T_obs', 'zscore', 'pvalue',
            '-log10pvalue', '-log10pvalue_smoothed',
        ])
    else:
        plot_df = pd.concat(plot_rows, ignore_index=True)

    # Merge in coverage / fold-change columns from the aggregation step so
    # downstream code that inspects the table keeps working.
    plot_df = plot_df.merge(
        agg[['name', 'position', 'coverage', 'number_exons', 'norm_coverage',
             'control_coverage', 'control_number_exons',
             'control_norm_coverage', 'fold_change']],
        on=['name', 'position'], how='left'
    )

    # Cluster detection / correction has been removed (results were not
    # actionable in practice); return an empty clusters table for backward
    # compatibility with the tuple-returning signature.
    clusters_df = pd.DataFrame(columns=[
        'name', 'label', 'start_pos', 'end_pos',
        'mass', 'cluster_pvalue', 'sign',
    ])

    return plot_df, clusters_df

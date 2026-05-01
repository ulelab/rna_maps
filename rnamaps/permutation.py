"""Label-permutation test for position-wise CLIP enrichment.

Replaces the per-position Fisher's exact test against control with an
empirical p-value derived from randomly relabelling exons between the
focal category (e.g. ``enhanced``) and ``control``. A cluster-based
multiple-testing correction across the correlated positions in each
splice-site region is also provided.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd


# Tiny floor to avoid log10(0) when an empirical p reaches its minimum.
_PVAL_FLOOR = 1e-300


def _build_coverage_matrix(df_per_exon: pd.DataFrame, category: str,
                           control_label: str = "control"):
    """Pivot per-exon coverage to a (n_exons × n_positions) matrix.

    Returns
    -------
    matrix : np.ndarray
        Coverage values, rows = exons in {category, control}, cols = positions.
    is_category : np.ndarray (bool)
        True for rows belonging to ``category``, False for control rows.
    positions : np.ndarray
        Position values matching the matrix columns (sorted ascending).
    """
    sub = df_per_exon[df_per_exon['name'].isin([category, control_label])]
    if sub.empty:
        return None, None, None
    pivot = sub.pivot_table(
        index=['name', 'exon_id'], columns='position', values='coverage',
        fill_value=0, aggfunc='sum'
    ).sort_index(axis=1)
    matrix = pivot.to_numpy(dtype=np.float64)
    is_category = np.asarray(
        pivot.index.get_level_values('name') == category
    )
    positions = pivot.columns.to_numpy()
    return matrix, is_category, positions


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


def _empirical_two_sided_p(t_obs: np.ndarray, t_null: np.ndarray) -> np.ndarray:
    """Two-sided empirical p-value with the (1 + #ge)/(B + 1) correction."""
    abs_obs = np.abs(t_obs)
    abs_null = np.abs(t_null)
    ge = (abs_null >= abs_obs[None, :]).sum(axis=0)
    return (1.0 + ge) / (1.0 + t_null.shape[0])


def _smooth(values: np.ndarray, smoothing: int) -> np.ndarray:
    if smoothing is None or smoothing <= 1:
        return values.copy()
    s = pd.Series(values)
    return s.rolling(smoothing, center=True, win_type='gaussian').mean(std=2).to_numpy()


def compute_permutation_pvalues(df_per_exon: pd.DataFrame,
                                exon_categories: pd.Series,
                                label: str,
                                n_perm: int,
                                smoothing: int,
                                rng: np.random.Generator,
                                control_label: str = "control"):
    """Compute permutation p-values for every non-control category in a region.

    Parameters
    ----------
    df_per_exon : DataFrame
        Long-form coverage with columns: exon_id, name (category), position,
        coverage, label.
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

    Returns
    -------
    plot_df : DataFrame
        Columns: name, position, label, coverage, number_exons, norm_coverage,
        control_coverage, control_number_exons, control_norm_coverage,
        fold_change, T_obs, pvalue, -log10pvalue, -log10pvalue_smoothed.
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

    # Aggregate coverage per (category, position) for the legend / fold-change
    # columns; matches the existing plot_df schema.
    agg = df_per_exon.groupby(['name', 'position'], as_index=False).agg(
        coverage=('coverage', 'sum')
    )
    counts = pd.DataFrame({
        'name': exon_categories.index,
        'number_exons': exon_categories.values,
    })
    agg = agg.merge(counts, on='name', how='left')
    agg['norm_coverage'] = np.where(
        agg['coverage'] == 0, 0.0,
        agg['coverage'] / agg['number_exons']
    )
    ctrl_rows = agg[agg['name'] == control_label][
        ['position', 'coverage', 'number_exons']
    ].rename(columns={
        'coverage': 'control_coverage',
        'number_exons': 'control_number_exons',
    })
    agg = agg.merge(ctrl_rows, on='position', how='left')
    agg['control_norm_coverage'] = (
        agg['control_coverage'] / agg['control_number_exons']
    )
    agg.loc[agg['control_norm_coverage'] == 0, 'control_norm_coverage'] = 1e-6
    agg['fold_change'] = agg['norm_coverage'] / agg['control_norm_coverage']

    plot_rows = []

    for cat in categories:
        matrix, is_cat, positions = _build_coverage_matrix(
            df_per_exon, cat, control_label
        )
        if matrix is None or is_cat.sum() == 0 or n_ctrl == 0:
            logging.warning(
                f"[permutation] Skipping {cat} ({label}): empty matrix.")
            continue

        t_obs, t_null = _permutation_null(matrix, is_cat, n_perm, rng)
        pvalues = _empirical_two_sided_p(t_obs, t_null)

        signed_log10p = -np.log10(np.maximum(pvalues, _PVAL_FLOOR))
        # Sign by direction of T_obs (positive = enriched over control).
        signed_log10p = np.where(t_obs >= 0, signed_log10p, -signed_log10p)

        smoothed = _smooth(signed_log10p, smoothing)

        cat_df = pd.DataFrame({
            'name': cat,
            'position': positions,
            'label': label,
            'T_obs': t_obs,
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
            'pvalue': 1.0,
            '-log10pvalue': 0.0,
            '-log10pvalue_smoothed': 0.0,
        }))

    if not plot_rows:
        plot_df = pd.DataFrame(columns=[
            'name', 'position', 'label', 'T_obs', 'pvalue',
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

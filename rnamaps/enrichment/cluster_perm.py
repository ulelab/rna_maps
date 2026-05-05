"""Cluster-based permutation test (Maris-Oostenveld 2007).

Per-region procedure for each non-control category vs. control:

1. Compute per-position Welch t-statistic ``t_obs(p)``.
2. Find contiguous runs where ``|t_obs(p)| > cluster_thresh``; cluster
   mass is ``sum(t_obs)`` over the run.
3. Permutation loop ``B = n_perm`` times: shuffle category/control
   labels, recompute Welch t, find clusters, record the maximum
   ``|cluster_mass|`` under that shuffle.
4. Cluster p-value:
   ``p_cluster = (1 + #(null_max >= |obs_mass|)) / (B + 1)``.

This controls family-wise error across positions while respecting the
position-to-position correlation that pointwise tests ignore. Unlike
the empirical pointwise p, cluster p still saturates at ``1/(B+1)``,
but in practice that ceiling is rarely hit because cluster mass is the
sum of an entire peak's t-values.
"""

from __future__ import annotations

import logging
from typing import List, Tuple

import numpy as np
import pandas as pd

from rnamaps.enrichment._base import EnrichmentResult
from rnamaps.permutation import _build_coverage_matrix


def _welch_t(cat_mat: np.ndarray, ctrl_mat: np.ndarray) -> np.ndarray:
    """Welch t-statistic per column, robust to zero-variance positions.

    Parameters
    ----------
    cat_mat : ndarray, shape (n_cat, n_pos)
    ctrl_mat : ndarray, shape (n_ctrl, n_pos)

    Returns
    -------
    ndarray, shape (n_pos,)
    """
    n_cat = cat_mat.shape[0]
    n_ctrl = ctrl_mat.shape[0]
    if n_cat < 2 or n_ctrl < 2:
        return np.zeros(cat_mat.shape[1], dtype=np.float64)
    mean_diff = cat_mat.mean(axis=0) - ctrl_mat.mean(axis=0)
    var_cat = cat_mat.var(axis=0, ddof=1)
    var_ctrl = ctrl_mat.var(axis=0, ddof=1)
    se = np.sqrt(var_cat / n_cat + var_ctrl / n_ctrl)
    t = np.zeros_like(mean_diff)
    valid = se > 0
    t[valid] = mean_diff[valid] / se[valid]
    return t


def _find_clusters(t: np.ndarray, thresh: float) -> List[Tuple[int, int, float, int]]:
    """Find contiguous suprathreshold runs in ``t``.

    Returns
    -------
    list of (start_idx, end_idx_inclusive, mass, sign)
        Sign is +1 if cluster mean is positive, -1 if negative.
        Both signs are returned (positive and negative clusters
        considered separately) so that a flip from + to - breaks the
        cluster.
    """
    if t.size == 0:
        return []
    sign = np.sign(t)
    above = np.abs(t) > thresh
    # Break clusters at sign flips so that a +/- transition starts a new run.
    n = t.size
    clusters: List[Tuple[int, int, float, int]] = []
    i = 0
    while i < n:
        if not above[i]:
            i += 1
            continue
        s = sign[i]
        j = i
        while j + 1 < n and above[j + 1] and sign[j + 1] == s:
            j += 1
        mass = float(t[i:j + 1].sum())
        clusters.append((i, j, mass, int(s) if s != 0 else 1))
        i = j + 1
    return clusters


def _max_abs_cluster_mass(t: np.ndarray, thresh: float) -> float:
    clusters = _find_clusters(t, thresh)
    if not clusters:
        return 0.0
    return max(abs(c[2]) for c in clusters)


def compute(
    df_per_exon: pd.DataFrame,
    exon_categories: pd.Series,
    label: str,
    *,
    rng: np.random.Generator,
    n_perm: int = 1000,
    cluster_thresh: float = 2.0,
    control_label: str = "control",
) -> EnrichmentResult:
    """Run cluster-based permutation test for each non-control category.

    Parameters
    ----------
    df_per_exon : DataFrame
        Long-form per-exon coverage with columns ``exon_id``, ``name``,
        ``position``, ``coverage``, ``label``.
    exon_categories : Series
        Counts per category.
    label : str
        Splice-site region label.
    rng : np.random.Generator
    n_perm : int
        Number of label permutations.
    cluster_thresh : float
        Cluster-defining ``|t|`` threshold.

    Returns
    -------
    EnrichmentResult
        ``plot_df`` columns: ``name, position, label, t_obs,
        in_sig_cluster, cluster_id, cluster_pvalue``.
        ``clusters_df`` columns: ``name, label, start_pos, end_pos, mass,
        cluster_pvalue, sign``.
    """
    if control_label not in exon_categories.index:
        raise ValueError(
            f"cluster_perm requires a '{control_label}' category."
        )

    categories = [c for c in exon_categories.index if c != control_label]
    n_ctrl_total = int(exon_categories.loc[control_label])

    plot_rows = []
    cluster_rows = []

    for cat in categories:
        matrix, is_cat, positions = _build_coverage_matrix(
            df_per_exon, cat, control_label
        )
        if matrix is None or is_cat.sum() == 0 or n_ctrl_total == 0:
            logging.warning(
                f"[cluster_perm] Skipping {cat} ({label}): empty matrix."
            )
            continue

        n_total = matrix.shape[0]
        n_c = int(is_cat.sum())

        cat_mat = matrix[is_cat]
        ctrl_mat = matrix[~is_cat]
        t_obs = _welch_t(cat_mat, ctrl_mat)
        obs_clusters = _find_clusters(t_obs, cluster_thresh)

        null_max = np.empty(n_perm, dtype=np.float64)
        for b in range(n_perm):
            idx = rng.choice(n_total, size=n_c, replace=False)
            mask = np.zeros(n_total, dtype=bool)
            mask[idx] = True
            t_b = _welch_t(matrix[mask], matrix[~mask])
            null_max[b] = _max_abs_cluster_mass(t_b, cluster_thresh)

        n_pos = positions.size
        in_sig = np.zeros(n_pos, dtype=bool)
        cluster_id = np.full(n_pos, -1, dtype=int)
        cluster_pvalue = np.ones(n_pos, dtype=np.float64)

        for k, (s, e, mass, sign) in enumerate(obs_clusters):
            p_clu = (1 + int((null_max >= abs(mass)).sum())) / (n_perm + 1)
            cluster_rows.append({
                'name': cat,
                'label': label,
                'start_pos': int(positions[s]),
                'end_pos': int(positions[e]),
                'mass': float(mass),
                'cluster_pvalue': float(p_clu),
                'sign': int(sign),
            })
            cluster_id[s:e + 1] = k
            cluster_pvalue[s:e + 1] = p_clu
            if p_clu <= 0.05:
                in_sig[s:e + 1] = True

        plot_rows.append(pd.DataFrame({
            'name': cat,
            'position': positions,
            'label': label,
            't_obs': t_obs,
            'in_sig_cluster': in_sig,
            'cluster_id': cluster_id,
            'cluster_pvalue': cluster_pvalue,
        }))

    # Control rows for legend continuity (zero stat).
    ctrl_positions = sorted(
        df_per_exon.loc[df_per_exon['name'] == control_label, 'position']
        .unique()
    )
    if ctrl_positions:
        plot_rows.append(pd.DataFrame({
            'name': control_label,
            'position': ctrl_positions,
            'label': label,
            't_obs': 0.0,
            'in_sig_cluster': False,
            'cluster_id': -1,
            'cluster_pvalue': 1.0,
        }))

    if not plot_rows:
        plot_df = pd.DataFrame(columns=[
            'name', 'position', 'label', 't_obs',
            'in_sig_cluster', 'cluster_id', 'cluster_pvalue',
        ])
    else:
        plot_df = pd.concat(plot_rows, ignore_index=True)

    clusters_df = pd.DataFrame(cluster_rows, columns=[
        'name', 'label', 'start_pos', 'end_pos',
        'mass', 'cluster_pvalue', 'sign',
    ])

    return EnrichmentResult(
        plot_df=plot_df,
        clusters_df=clusters_df,
        plot_kind='clusters',
        y_columns=('t_obs',),
        method_name='cluster_perm',
        ylabel='Welch t (cluster-permutation)',
    )

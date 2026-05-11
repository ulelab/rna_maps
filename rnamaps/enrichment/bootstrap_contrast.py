"""Bootstrap-based contrast estimators for category vs. control.

For each non-control category vs. control, exons are resampled with
replacement ``B`` times and two contrasts are reported per position with
their 2.5/97.5 percentile bands:

- ``delta(p)  = mean_cov_cat(p) - mean_cov_ctrl(p)`` (additive scale).
- ``log2fc(p) = log2((mean_cov_cat(p) + eps) / (mean_cov_ctrl(p) + eps))``
  (multiplicative scale, library-size invariant).

The bootstrap reflects sampling uncertainty in the estimator. It does
*not* correct for bias from contaminated controls (silently regulated
exons in the control pool). Use the ``--control_set`` flag for that.

When ``smoothing > 1``, each bootstrap iteration's per-position contrast
is convolved with a centred Gaussian-weighted rolling mean (matching the
kernel in ``rnamaps.permutation._smooth``) *before* taking the across-
iteration mean and percentiles. Smoothing first makes the CI band the
correct uncertainty band for the smoothed estimator; smoothing the
percentiles afterwards would understate uncertainty at sharp features.
"""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd

from rnamaps.enrichment._base import EnrichmentResult
from rnamaps.permutation import _build_coverage_matrix


def _smooth_along_positions(
    arr: np.ndarray, smoothing: int
) -> np.ndarray:
    """Centred Gaussian rolling-mean smoothing along the position axis.

    Matches ``rnamaps.permutation._smooth`` (window ``smoothing``, Gaussian
    weights with ``std=2``). Accepts a 1D ``(n_pos,)`` array or a 2D
    ``(n_boot, n_pos)`` array; for the 2D case each bootstrap row is
    smoothed independently in vectorised C via ``DataFrame.rolling``.

    Edge positions where the window does not fully fit become NaN, mirroring
    the existing permutation-test behaviour.
    """
    if smoothing is None or smoothing <= 1:
        return arr
    if arr.ndim == 1:
        s = pd.Series(arr)
        return s.rolling(
            smoothing, center=True, win_type='gaussian'
        ).mean(std=2).to_numpy()
    # 2D: transpose so positions are rows and each bootstrap iteration is a
    # column; pandas rolls along axis 0 (the default) per-column → smooths
    # each bootstrap row in vectorised C. Transpose back at the end.
    df = pd.DataFrame(arr.T)
    smoothed = df.rolling(
        smoothing, center=True, win_type='gaussian'
    ).mean(std=2).to_numpy()
    return smoothed.T


def _bootstrap_means(
    matrix: np.ndarray,
    n: int,
    n_boot: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Bootstrap mean coverage per position from `matrix`.

    Parameters
    ----------
    matrix : ndarray, shape (n_rows, n_positions)
        Per-exon coverage rows.
    n : int
        Number of rows to draw with replacement per iteration (usually
        ``matrix.shape[0]``).
    n_boot : int
    rng : np.random.Generator

    Returns
    -------
    ndarray, shape (n_boot, n_positions)
    """
    if matrix.ndim != 2 or matrix.shape[0] == 0:
        n_pos = matrix.shape[1] if matrix.ndim == 2 else 0
        return np.zeros((n_boot, n_pos), dtype=np.float64)
    out = np.empty((n_boot, matrix.shape[1]), dtype=np.float64)
    n_rows = matrix.shape[0]
    for b in range(n_boot):
        idx = rng.integers(0, n_rows, size=n)
        out[b] = matrix[idx].mean(axis=0)
    return out


def _aggregate_legacy_columns(
    df_per_exon: pd.DataFrame,
    exon_categories: pd.Series,
    control_label: str,
) -> pd.DataFrame:
    """Recreate the coverage / fold_change columns the legacy schema used.

    Returns a DataFrame indexed implicitly by (name, position) with the
    same columns the existing pipeline TSV exposes for backward compat.
    """
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
    return agg


def compute(
    df_per_exon: pd.DataFrame,
    exon_categories: pd.Series,
    label: str,
    *,
    rng: np.random.Generator,
    n_boot: int = 1000,
    pseudocount: Optional[float] = None,
    pseudocount_frac: float = 0.01,
    bootstrap_control_fixed: bool = False,
    smoothing: int = 1,
    control_label: str = "control",
    ci_low: float = 2.5,
    ci_high: float = 97.5,
) -> EnrichmentResult:
    """Compute bootstrap-CI contrasts (delta and log2fc) for each non-control category.

    Parameters
    ----------
    df_per_exon : DataFrame
        Long-form per-exon coverage with columns: ``exon_id``, ``name``,
        ``position``, ``coverage``, ``label``.
    exon_categories : Series
        Counts per category.
    label : str
        Splice-site region label, e.g. ``"middle_3ss"``.
    rng : np.random.Generator
    n_boot : int
        Bootstrap iterations.
    pseudocount : float or None
        If not None, fixed pseudocount for log2fc. Otherwise an adaptive
        ``eps = max(1e-3, pseudocount_frac * median(mean_cov_ctrl over region))``
        is used per region/category to keep log2fc stable at sparse positions.
    pseudocount_frac : float
        Fraction of regional median ctrl coverage when adaptive.
    bootstrap_control_fixed : bool
        If True, treat the control mean as a constant and skip resampling
        the control group. Equivalent up to negligible variance when
        ``n_ctrl >= 2000`` and ~5-10x faster.
    smoothing : int
        Centred Gaussian rolling-mean window (in positions) applied to
        each bootstrap iteration's ``delta`` / ``log2fc`` before taking
        the across-iteration mean and percentiles. ``smoothing <= 1``
        disables smoothing. Matches the kernel used by ``permutation_z``.
    control_label : str
    ci_low, ci_high : float
        Percentile bounds for the CI band.

    Returns
    -------
    EnrichmentResult
        ``plot_df`` carries ``delta``, ``delta_lo``, ``delta_hi``,
        ``log2fc``, ``log2fc_lo``, ``log2fc_hi`` plus the legacy
        coverage / fold_change parity columns.
    """
    if control_label not in exon_categories.index:
        raise ValueError(
            f"bootstrap_contrast requires a '{control_label}' category."
        )

    categories = [c for c in exon_categories.index if c != control_label]
    n_ctrl_total = int(exon_categories.loc[control_label])

    agg = _aggregate_legacy_columns(df_per_exon, exon_categories, control_label)

    plot_rows = []

    for cat in categories:
        matrix, is_cat, positions = _build_coverage_matrix(
            df_per_exon, cat, control_label
        )
        if matrix is None or is_cat.sum() == 0 or n_ctrl_total == 0:
            logging.warning(
                f"[bootstrap_contrast] Skipping {cat} ({label}): empty matrix."
            )
            continue

        cat_mat = matrix[is_cat]
        ctrl_mat = matrix[~is_cat]
        n_c = cat_mat.shape[0]
        n_ctrl = ctrl_mat.shape[0]

        if not bootstrap_control_fixed and n_c > 0 and n_ctrl >= 20 * n_c:
            logging.info(
                f"[bootstrap_contrast] {cat} ({label}): n_ctrl={n_ctrl} "
                f">= 20x n_{cat}={n_c}. --bootstrap_control_fixed would be "
                f"~{n_ctrl / n_c:.0f}x faster with negligible CI change."
            )

        ctrl_mean = ctrl_mat.mean(axis=0) if n_ctrl > 0 else np.zeros(
            matrix.shape[1]
        )
        if pseudocount is not None:
            eps_ = float(pseudocount)
        else:
            med = float(np.median(ctrl_mean)) if n_ctrl > 0 else 0.0
            eps_ = max(1e-3, pseudocount_frac * med)

        cat_means = _bootstrap_means(cat_mat, n_c, n_boot, rng)
        if bootstrap_control_fixed:
            ctrl_means = np.broadcast_to(
                ctrl_mean, (n_boot, ctrl_mean.size)
            )
        else:
            ctrl_means = _bootstrap_means(ctrl_mat, n_ctrl, n_boot, rng)

        delta_b = cat_means - ctrl_means
        log2fc_b = np.log2((cat_means + eps_) / (ctrl_means + eps_))

        # Smooth each bootstrap iteration along the position axis before
        # collapsing to mean/percentiles so the CI band is the uncertainty
        # of the smoothed estimator, not the smoothed uncertainty.
        delta_b = _smooth_along_positions(delta_b, smoothing)
        log2fc_b = _smooth_along_positions(log2fc_b, smoothing)

        cat_df = pd.DataFrame({
            'name': cat,
            'position': positions,
            'label': label,
            'delta': np.nanmean(delta_b, axis=0),
            'delta_lo': np.nanpercentile(delta_b, ci_low, axis=0),
            'delta_hi': np.nanpercentile(delta_b, ci_high, axis=0),
            'log2fc': np.nanmean(log2fc_b, axis=0),
            'log2fc_lo': np.nanpercentile(log2fc_b, ci_low, axis=0),
            'log2fc_hi': np.nanpercentile(log2fc_b, ci_high, axis=0),
            'pseudocount': eps_,
        })
        plot_rows.append(cat_df)

    ctrl_positions = agg.loc[agg['name'] == control_label, 'position'].to_numpy()
    if len(ctrl_positions):
        plot_rows.append(pd.DataFrame({
            'name': control_label,
            'position': ctrl_positions,
            'label': label,
            'delta': 0.0,
            'delta_lo': 0.0,
            'delta_hi': 0.0,
            'log2fc': 0.0,
            'log2fc_lo': 0.0,
            'log2fc_hi': 0.0,
            'pseudocount': np.nan,
        }))

    if not plot_rows:
        plot_df = pd.DataFrame(columns=[
            'name', 'position', 'label',
            'delta', 'delta_lo', 'delta_hi',
            'log2fc', 'log2fc_lo', 'log2fc_hi',
            'pseudocount',
        ])
    else:
        plot_df = pd.concat(plot_rows, ignore_index=True)

    plot_df = plot_df.merge(
        agg[['name', 'position', 'coverage', 'number_exons', 'norm_coverage',
             'control_coverage', 'control_number_exons',
             'control_norm_coverage', 'fold_change']],
        on=['name', 'position'], how='left'
    )

    return EnrichmentResult(
        plot_df=plot_df,
        clusters_df=pd.DataFrame(),
        plot_kind='ribbon',
        y_columns=('delta', 'log2fc'),
        method_name='bootstrap_contrast',
        ylabel='bootstrap contrast vs control',
    )

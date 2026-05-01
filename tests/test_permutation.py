"""Tests for the label-permutation test in rnamaps.permutation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rnamaps.permutation import (
    _empirical_two_sided_p,
    _permutation_null,
    compute_permutation_pvalues,
)


def _make_df(matrix, names, positions, label="middle_3ss"):
    """Build the long-form per-exon coverage frame the tests need."""
    rows = []
    for i, (name, row) in enumerate(zip(names, matrix)):
        for pos, cov in zip(positions, row):
            rows.append({
                'exon_id': f'e{i}',
                'name': name,
                'position': int(pos),
                'coverage': float(cov),
                'label': label,
            })
    return pd.DataFrame(rows)


def test_null_type_i_error_controlled():
    """Under H0 the empirical false-positive rate at alpha=0.05 is <= ~0.05.

    The permutation p-value is exact and conservative under exchangeability
    (especially for discrete data with ties), so the rejection rate must
    not exceed the nominal alpha by more than sampling error.
    """
    rng = np.random.default_rng(0)
    n_exons, n_pos = 60, 40
    matrix = rng.poisson(lam=0.4, size=(n_exons, n_pos)).astype(float)
    is_cat = np.zeros(n_exons, dtype=bool)
    is_cat[:30] = True

    _, t_null = _permutation_null(matrix, is_cat, n_perm=200, rng=rng)
    # Treat each permutation in turn as the "observed" and test against the
    # rest -- a standard self-consistency check for permutation calibration.
    rejections = 0
    total = 0
    for b in range(t_null.shape[0]):
        t_obs = t_null[b]
        t_rest = np.delete(t_null, b, axis=0)
        pvals = _empirical_two_sided_p(t_obs, t_rest)
        rejections += int((pvals <= 0.05).sum())
        total += pvals.size
    rate = rejections / total
    # Allow some slack: 0.05 nominal, with ~8000 tests SE ~ 0.0024.
    assert rate <= 0.06, f"Type-I error rate too high: {rate:.4f}"


def test_strong_signal_recovered():
    """A position with large category-vs-control gap gets the smallest p."""
    rng = np.random.default_rng(1)
    n_pos = 20
    spike_pos = 7
    cat = np.zeros((30, n_pos))
    ctrl = np.zeros((30, n_pos))
    cat[:, spike_pos] = 5.0  # massive enrichment at one position
    matrix = np.vstack([cat, ctrl])
    is_cat = np.array([True] * 30 + [False] * 30)

    t_obs, t_null = _permutation_null(matrix, is_cat, n_perm=500, rng=rng)
    pvals = _empirical_two_sided_p(t_obs, t_null)
    assert int(np.argmin(pvals)) == spike_pos
    assert pvals[spike_pos] == pytest.approx(1.0 / 501, rel=1e-6)


def test_seed_determinism():
    """Same seed -> identical p-values; different seed -> different."""
    n_exons, n_pos = 40, 25
    rng_data = np.random.default_rng(7)
    matrix = rng_data.poisson(lam=0.6, size=(n_exons, n_pos)).astype(float)
    is_cat = np.zeros(n_exons, dtype=bool)
    is_cat[:20] = True

    rng_a = np.random.default_rng(123)
    rng_b = np.random.default_rng(123)
    rng_c = np.random.default_rng(124)
    _, null_a = _permutation_null(matrix, is_cat, n_perm=100, rng=rng_a)
    _, null_b = _permutation_null(matrix, is_cat, n_perm=100, rng=rng_b)
    _, null_c = _permutation_null(matrix, is_cat, n_perm=100, rng=rng_c)
    np.testing.assert_array_equal(null_a, null_b)
    assert not np.array_equal(null_a, null_c)


def test_compute_permutation_pvalues_end_to_end():
    """Full DataFrame in/out works and returns the expected schema."""
    rng = np.random.default_rng(2)
    n_pos = 15
    spike_pos = 5
    cat = np.zeros((20, n_pos))
    ctrl = np.zeros((20, n_pos))
    cat[:, spike_pos] = 3.0
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * 20 + ['control'] * 20
    positions = np.arange(n_pos)
    df = _make_df(matrix, names, positions, label='middle_3ss')

    exon_categories = pd.Series({'enhanced': 20, 'control': 20})
    plot_df, clusters_df = compute_permutation_pvalues(
        df, exon_categories, label='middle_3ss',
        n_perm=200, smoothing=1, rng=rng,
    )
    # Schema
    expected_cols = {
        'name', 'position', 'label', 'T_obs', 'pvalue',
        '-log10pvalue', '-log10pvalue_smoothed',
        'coverage', 'number_exons', 'norm_coverage',
        'control_coverage', 'control_number_exons',
        'control_norm_coverage', 'fold_change',
    }
    assert expected_cols.issubset(plot_df.columns)
    # Spike position must be the most significant for "enhanced"
    enh = plot_df[plot_df['name'] == 'enhanced'].set_index('position')
    assert enh['pvalue'].idxmin() == spike_pos
    # clusters_df is now always empty (cluster output removed).
    assert clusters_df.empty

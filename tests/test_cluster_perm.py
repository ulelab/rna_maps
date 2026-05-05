"""Tests for rnamaps.enrichment.cluster_perm."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rnamaps.enrichment.cluster_perm import (
    _find_clusters,
    _welch_t,
    compute,
)


def _make_df(matrix, names, positions, label="middle_3ss"):
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


def test_find_clusters_basic():
    t = np.array([0.0, 3.0, 4.0, 0.5, -3.5, -2.5, 0.0, 5.0])
    clusters = _find_clusters(t, thresh=2.0)
    starts = [c[0] for c in clusters]
    ends = [c[1] for c in clusters]
    masses = [c[2] for c in clusters]
    signs = [c[3] for c in clusters]
    assert starts == [1, 4, 7]
    assert ends == [2, 5, 7]
    assert signs == [1, -1, 1]
    np.testing.assert_allclose(masses, [7.0, -6.0, 5.0])


def test_find_clusters_sign_break():
    """Adjacent +/- positions don't merge into one cluster."""
    t = np.array([3.0, -3.0])
    clusters = _find_clusters(t, thresh=2.0)
    assert len(clusters) == 2


def test_welch_t_zero_when_equal():
    cat = np.ones((5, 4))
    ctrl = np.ones((5, 4))
    t = _welch_t(cat, ctrl)
    np.testing.assert_array_equal(t, np.zeros(4))


def test_recovers_planted_peak():
    """A wide Gaussian peak in cat (vs flat ctrl) yields a sig cluster."""
    rng = np.random.default_rng(0)
    n_pos = 60
    n_cat = 40
    n_ctrl = 40
    centre = 30
    width = 6
    x = np.arange(n_pos)
    bump = 2.0 * np.exp(-0.5 * ((x - centre) / width) ** 2)
    cat = rng.normal(scale=0.5, size=(n_cat, n_pos)) + bump
    ctrl = rng.normal(scale=0.5, size=(n_ctrl, n_pos))
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * n_cat + ['control'] * n_ctrl
    df = _make_df(matrix, names, np.arange(n_pos))

    exon_categories = pd.Series({'enhanced': n_cat, 'control': n_ctrl})
    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_perm=200, cluster_thresh=2.0)
    cl = res.clusters_df[res.clusters_df['name'] == 'enhanced']
    assert not cl.empty, "No clusters detected for planted peak"
    # The largest cluster should overlap the peak centre.
    big = cl.iloc[cl['mass'].abs().argmax()]
    assert big['start_pos'] <= centre <= big['end_pos']
    # And it should be highly significant.
    assert big['cluster_pvalue'] < 0.05
    assert big['sign'] == 1


def test_null_data_few_significant_clusters():
    """Under H0, the rate of any-significant cluster should stay low."""
    rng = np.random.default_rng(2)
    n_pos = 40
    n_cat = 30
    n_ctrl = 30
    matrix = rng.normal(size=(n_cat + n_ctrl, n_pos))
    names = ['enhanced'] * n_cat + ['control'] * n_ctrl
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': n_cat, 'control': n_ctrl})

    n_trials = 30
    sig = 0
    for trial in range(n_trials):
        rng_t = np.random.default_rng(1000 + trial)
        # Re-sample data per trial to vary the null instances.
        m_t = rng_t.normal(size=(n_cat + n_ctrl, n_pos))
        df_t = _make_df(m_t, names, np.arange(n_pos))
        res = compute(df_t, exon_categories, label='middle_3ss',
                      rng=rng_t, n_perm=200, cluster_thresh=2.0)
        cl = res.clusters_df[res.clusters_df['name'] == 'enhanced']
        if (cl['cluster_pvalue'] <= 0.05).any():
            sig += 1
    rate = sig / n_trials
    assert rate <= 0.20, f"FWER under H0 too high: {rate:.2f}"


def test_result_schema_and_metadata():
    rng = np.random.default_rng(3)
    matrix = rng.normal(size=(20, 8))
    names = ['enhanced'] * 10 + ['control'] * 10
    df = _make_df(matrix, names, np.arange(8))
    exon_categories = pd.Series({'enhanced': 10, 'control': 10})
    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_perm=100)
    expected = {'name', 'position', 'label', 't_obs',
                'in_sig_cluster', 'cluster_id', 'cluster_pvalue'}
    assert expected.issubset(res.plot_df.columns)
    assert res.plot_kind == 'clusters'
    assert res.method_name == 'cluster_perm'


def test_missing_control_raises():
    rng = np.random.default_rng(7)
    df = _make_df(np.zeros((4, 3)), ['enhanced'] * 4, np.arange(3))
    exon_categories = pd.Series({'enhanced': 4})
    with pytest.raises(ValueError, match="control"):
        compute(df, exon_categories, label='middle_3ss',
                rng=rng, n_perm=10)

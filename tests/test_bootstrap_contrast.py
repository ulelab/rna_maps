"""Tests for rnamaps.enrichment.bootstrap_contrast."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rnamaps.enrichment.bootstrap_contrast import (
    _bootstrap_means,
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


def test_bootstrap_means_recovers_mean():
    rng = np.random.default_rng(0)
    matrix = rng.normal(size=(200, 10))
    boot = _bootstrap_means(matrix, n=200, n_boot=500, rng=rng)
    # Bootstrap mean of mean estimates should track the sample mean.
    np.testing.assert_allclose(
        boot.mean(axis=0), matrix.mean(axis=0), atol=0.05
    )


def test_strong_signal_recovered():
    """A spike position should give a positive log2fc with CI excluding 0."""
    rng = np.random.default_rng(1)
    n_pos = 12
    spike = 4
    cat = np.zeros((30, n_pos))
    ctrl = np.zeros((30, n_pos))
    cat[:, spike] = 5.0
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * 30 + ['control'] * 30
    positions = np.arange(n_pos)
    df = _make_df(matrix, names, positions)

    exon_categories = pd.Series({'enhanced': 30, 'control': 30})
    res = compute(
        df, exon_categories, label='middle_3ss',
        rng=rng, n_boot=500,
    )
    out = res.plot_df
    enh = out[out['name'] == 'enhanced'].set_index('position')

    # log2fc and delta must be very positive at the spike.
    assert enh.loc[spike, 'log2fc'] > 5
    assert enh.loc[spike, 'delta'] > 4
    # CI lower bound at spike must exclude zero.
    assert enh.loc[spike, 'log2fc_lo'] > 0
    assert enh.loc[spike, 'delta_lo'] > 0
    # Off-spike positions: contrast hugs zero.
    other = enh.drop(spike)
    assert other['delta'].abs().max() < 0.5


def test_null_data_ci_brackets_zero():
    """Under H0 the CI for delta should generally include 0."""
    rng = np.random.default_rng(2)
    n_exons, n_pos = 60, 15
    matrix = rng.poisson(lam=0.5, size=(n_exons, n_pos)).astype(float)
    names = ['enhanced'] * 30 + ['control'] * 30
    positions = np.arange(n_pos)
    df = _make_df(matrix, names, positions)

    exon_categories = pd.Series({'enhanced': 30, 'control': 30})
    res = compute(
        df, exon_categories, label='middle_3ss',
        rng=rng, n_boot=500,
    )
    out = res.plot_df
    enh = out[out['name'] == 'enhanced']
    brackets = ((enh['delta_lo'] <= 0) & (enh['delta_hi'] >= 0)).mean()
    assert brackets >= 0.85, f"Only {brackets:.2f} of CIs bracket 0 under H0"


def test_class_imbalance_fixed_control_matches_full():
    """With n_ctrl >> n_cat, --bootstrap_control_fixed ~= full bootstrap."""
    rng_a = np.random.default_rng(3)
    rng_b = np.random.default_rng(3)
    n_pos = 8
    n_cat = 50
    n_ctrl = 2000
    rng_data = np.random.default_rng(99)
    cat = rng_data.poisson(lam=0.4, size=(n_cat, n_pos)).astype(float)
    ctrl = rng_data.poisson(lam=0.4, size=(n_ctrl, n_pos)).astype(float)
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * n_cat + ['control'] * n_ctrl
    positions = np.arange(n_pos)
    df = _make_df(matrix, names, positions)
    exon_categories = pd.Series({'enhanced': n_cat, 'control': n_ctrl})

    full = compute(df, exon_categories, label='middle_3ss',
                   rng=rng_a, n_boot=500).plot_df
    fixed = compute(df, exon_categories, label='middle_3ss',
                    rng=rng_b, n_boot=500,
                    bootstrap_control_fixed=True).plot_df

    full_e = full[full['name'] == 'enhanced'].set_index('position')
    fixed_e = fixed[fixed['name'] == 'enhanced'].set_index('position')

    # Point estimates should be close.
    np.testing.assert_allclose(
        fixed_e['delta'].values, full_e['delta'].values, atol=0.05
    )
    # CI widths should be within a small factor.
    full_w = (full_e['delta_hi'] - full_e['delta_lo']).values
    fixed_w = (fixed_e['delta_hi'] - fixed_e['delta_lo']).values
    ratio = fixed_w.mean() / full_w.mean()
    assert 0.6 <= ratio <= 1.4, f"CI width ratio {ratio:.2f} out of range"


def test_pseudocount_modes():
    """Adaptive pseudocount stays positive and override is respected."""
    rng = np.random.default_rng(4)
    n_pos = 6
    matrix = rng.poisson(lam=0.05, size=(40, n_pos)).astype(float)
    names = ['enhanced'] * 20 + ['control'] * 20
    positions = np.arange(n_pos)
    df = _make_df(matrix, names, positions)
    exon_categories = pd.Series({'enhanced': 20, 'control': 20})

    res_adaptive = compute(df, exon_categories, label='middle_3ss',
                           rng=rng, n_boot=200)
    enh = res_adaptive.plot_df[res_adaptive.plot_df['name'] == 'enhanced']
    assert (enh['pseudocount'] > 0).all()

    res_fixed = compute(df, exon_categories, label='middle_3ss',
                        rng=np.random.default_rng(4), n_boot=200,
                        pseudocount=0.05)
    enh2 = res_fixed.plot_df[res_fixed.plot_df['name'] == 'enhanced']
    assert np.allclose(enh2['pseudocount'].dropna().values, 0.05)


def test_result_schema_and_metadata():
    rng = np.random.default_rng(5)
    matrix = rng.poisson(lam=0.5, size=(20, 5)).astype(float)
    names = ['enhanced'] * 10 + ['control'] * 10
    positions = np.arange(5)
    df = _make_df(matrix, names, positions)
    exon_categories = pd.Series({'enhanced': 10, 'control': 10})
    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_boot=100)
    expected = {
        'name', 'position', 'label',
        'delta', 'delta_lo', 'delta_hi',
        'log2fc', 'log2fc_lo', 'log2fc_hi',
        'pseudocount', 'coverage', 'number_exons', 'norm_coverage',
        'control_coverage', 'control_number_exons',
        'control_norm_coverage', 'fold_change',
    }
    assert expected.issubset(res.plot_df.columns)
    assert res.plot_kind == 'ribbon'
    assert res.method_name == 'bootstrap_contrast'
    assert set(res.y_columns) == {'delta', 'log2fc'}


def test_smoothing_widens_no_signal_curve():
    """Heavy smoothing should reduce per-position fluctuations under H0."""
    rng = np.random.default_rng(8)
    matrix = rng.poisson(lam=0.5, size=(60, 30)).astype(float)
    names = ['enhanced'] * 30 + ['control'] * 30
    positions = np.arange(30)
    df = _make_df(matrix, names, positions)
    exon_categories = pd.Series({'enhanced': 30, 'control': 30})

    raw = compute(df, exon_categories, label='middle_3ss',
                  rng=np.random.default_rng(8), n_boot=200,
                  smoothing=1).plot_df
    smoothed = compute(df, exon_categories, label='middle_3ss',
                       rng=np.random.default_rng(8), n_boot=200,
                       smoothing=11).plot_df
    raw_e = raw[raw['name'] == 'enhanced']['delta'].to_numpy()
    sm_e = smoothed[smoothed['name'] == 'enhanced']['delta'].dropna().to_numpy()

    assert np.nanstd(sm_e) < np.nanstd(raw_e), (
        "Smoothed delta should fluctuate less than raw delta under H0."
    )


def test_missing_control_raises():
    rng = np.random.default_rng(6)
    df = _make_df(np.zeros((4, 3)), ['enhanced'] * 4, np.arange(3))
    exon_categories = pd.Series({'enhanced': 4})
    with pytest.raises(ValueError, match="control"):
        compute(df, exon_categories, label='middle_3ss',
                rng=rng, n_boot=10)

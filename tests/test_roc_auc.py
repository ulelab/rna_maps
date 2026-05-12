"""Tests for rnamaps.enrichment.roc_auc."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rnamaps.coverage import RegionCoverage
from rnamaps.enrichment.roc_auc import (
    _aggregate_window,
    _rank_auc,
    _roc_curve,
    _vectorised_auc_per_column,
    compute,
)


def _make_region(matrix, names, label="middle_3ss"):
    """Wrap a (n_exons, n_positions) matrix as a RegionCoverage."""
    n_pos = matrix.shape[1]
    return RegionCoverage(
        label=label,
        matrix=np.asarray(matrix, dtype=np.float64),
        exon_ids=np.array([f"e{i}" for i in range(len(names))],
                          dtype=object),
        exon_names=np.asarray(names, dtype=object),
        positions=np.arange(1, n_pos + 1, dtype=np.int32),
    )


def test_rank_auc_perfect_separation():
    """A perfect predictor (positives strictly above negatives) gives AUC=1."""
    scores = np.array([10.0, 9.0, 8.0, 1.0, 0.5, 0.0])
    labels = np.array([1, 1, 1, 0, 0, 0])
    assert _rank_auc(scores, labels) == pytest.approx(1.0)


def test_rank_auc_pure_noise():
    """Random labels at fixed scores give AUC ≈ 0.5 on average."""
    rng = np.random.default_rng(0)
    scores = rng.normal(size=2000)
    labels = rng.integers(0, 2, size=2000)
    auc = _rank_auc(scores, labels)
    assert 0.45 < auc < 0.55


def test_rank_auc_anti_correlated():
    """Inverted ranking gives AUC=0 (the symmetric counterpart)."""
    scores = np.array([0.1, 0.2, 0.3, 0.8, 0.9, 1.0])
    labels = np.array([1, 1, 1, 0, 0, 0])
    assert _rank_auc(scores, labels) == pytest.approx(0.0)


def test_rank_auc_degenerate_groups():
    """All-positive or all-negative inputs collapse to AUC=0.5."""
    scores = np.linspace(0, 1, 10)
    assert _rank_auc(scores, np.ones(10, dtype=int)) == 0.5
    assert _rank_auc(scores, np.zeros(10, dtype=int)) == 0.5


def test_vectorised_auc_matches_per_column_rank_auc():
    """Vectorised path agrees with the scalar Mann–Whitney identity."""
    rng = np.random.default_rng(1)
    n_exons, n_pos = 50, 30
    matrix = rng.normal(size=(n_exons, n_pos))
    is_pos = np.zeros(n_exons, dtype=bool)
    is_pos[:20] = True
    auc_vec = _vectorised_auc_per_column(matrix, is_pos)
    auc_loop = np.array(
        [_rank_auc(matrix[:, p], is_pos.astype(int))
         for p in range(n_pos)]
    )
    np.testing.assert_allclose(auc_vec, auc_loop, atol=1e-12)


def test_roc_curve_endpoints_and_monotonic():
    """ROC curve runs from (0,0) to (1,1) and is monotone non-decreasing."""
    rng = np.random.default_rng(2)
    n = 80
    labels = rng.integers(0, 2, size=n)
    scores = rng.normal(size=n) + labels * 0.7
    fpr, tpr = _roc_curve(scores, labels)
    assert fpr[0] == 0.0 and tpr[0] == 0.0
    assert fpr[-1] == 1.0 and tpr[-1] == 1.0
    assert np.all(np.diff(fpr) >= -1e-12)
    assert np.all(np.diff(tpr) >= -1e-12)


def test_aggregate_window_modes():
    matrix = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
    np.testing.assert_array_equal(
        _aggregate_window(matrix, "mean"), [1.0, 4.0]
    )
    np.testing.assert_array_equal(
        _aggregate_window(matrix, "max"), [2.0, 5.0]
    )
    np.testing.assert_array_equal(
        _aggregate_window(matrix, "sum"), [3.0, 12.0]
    )


def test_compute_strong_localised_signal():
    """A category-only spike at one position gives AUC ≈ 1 there and
    ≈ 0.5 elsewhere, and a region AUC well above 0.5 for both
    aggregators."""
    rng = np.random.default_rng(3)
    n_pos = 21
    spike = 10
    n_cat, n_ctrl = 40, 60
    matrix = np.zeros((n_cat + n_ctrl, n_pos))
    matrix[:n_cat, spike] = 1.0
    matrix += rng.normal(scale=0.01, size=matrix.shape)
    names = (["enhanced"] * n_cat) + (["control"] * n_ctrl)
    region = _make_region(matrix, names)
    exon_categories = pd.Series({"enhanced": n_cat, "control": n_ctrl})

    res = compute(
        region, exon_categories, "middle_3ss",
        rng=np.random.default_rng(7),
        roc_aggregator="both",
        n_perm=0,
        smoothing=1,
        binarise=False,
    )

    enh_pos = res.plot_df[res.plot_df['name'] == 'enhanced']
    enh_pos = enh_pos.sort_values('position')
    aucs = enh_pos['auc'].to_numpy()
    assert aucs[spike] > 0.95
    other = np.delete(aucs, spike)
    assert np.median(other) == pytest.approx(0.5, abs=0.05)

    region_aucs = res.extras['region_auc_df']
    enh_region = region_aucs[region_aucs['name'] == 'enhanced']
    assert (enh_region['auc'] > 0.7).all()

    # ROC curve is sane (starts at origin, ends at (1,1))
    roc = res.extras['roc_curves_df']
    enh_roc = roc[(roc['name'] == 'enhanced') & (roc['aggregator'] == 'mean')]
    assert enh_roc['fpr'].iloc[0] == 0.0
    assert enh_roc['tpr'].iloc[0] == 0.0
    assert enh_roc['fpr'].iloc[-1] == 1.0
    assert enh_roc['tpr'].iloc[-1] == 1.0


def test_compute_null_yields_auc_around_half():
    """When category and control are drawn from the same distribution
    the per-position AUC distribution centres on 0.5."""
    rng = np.random.default_rng(4)
    n_pos = 20
    n_cat, n_ctrl = 30, 30
    matrix = rng.normal(size=(n_cat + n_ctrl, n_pos))
    names = (["enhanced"] * n_cat) + (["control"] * n_ctrl)
    region = _make_region(matrix, names)
    exon_categories = pd.Series({"enhanced": n_cat, "control": n_ctrl})

    res = compute(
        region, exon_categories, "middle_3ss",
        rng=np.random.default_rng(7),
        roc_aggregator="mean",
        n_perm=0,
        smoothing=1,
        binarise=False,
    )
    enh = res.plot_df[res.plot_df['name'] == 'enhanced']
    assert np.median(enh['auc'].to_numpy()) == pytest.approx(0.5, abs=0.08)


def test_compute_permutation_pvalues_present():
    """With ``n_perm > 0`` per-position and per-region p-values are
    populated and lie in (0, 1]."""
    rng = np.random.default_rng(5)
    n_pos = 12
    n_cat, n_ctrl = 20, 20
    matrix = rng.normal(size=(n_cat + n_ctrl, n_pos))
    matrix[:n_cat, 5] += 1.5
    names = (["enhanced"] * n_cat) + (["control"] * n_ctrl)
    region = _make_region(matrix, names)
    exon_categories = pd.Series({"enhanced": n_cat, "control": n_ctrl})

    res = compute(
        region, exon_categories, "middle_3ss",
        rng=np.random.default_rng(7),
        roc_aggregator="mean",
        n_perm=200,
        smoothing=1,
        binarise=False,
    )
    enh = res.plot_df[res.plot_df['name'] == 'enhanced']
    pvals = enh['pvalue'].to_numpy()
    assert np.isfinite(pvals).all()
    assert ((pvals > 0) & (pvals <= 1)).all()

    region_p = res.extras['region_auc_df']
    assert np.isfinite(region_p['pvalue'].to_numpy()).all()


def test_compute_signed_auc_centred_on_zero_for_control():
    """Control rows should carry signed AUC = 0 (baseline)."""
    rng = np.random.default_rng(6)
    matrix = rng.normal(size=(40, 8))
    names = (["enhanced"] * 20) + (["control"] * 20)
    region = _make_region(matrix, names)
    exon_categories = pd.Series({"enhanced": 20, "control": 20})

    res = compute(
        region, exon_categories, "middle_3ss",
        rng=np.random.default_rng(7),
        roc_aggregator="mean", n_perm=0, smoothing=1, binarise=False,
    )
    ctrl = res.plot_df[res.plot_df['name'] == 'control']
    np.testing.assert_array_equal(ctrl['auc_signed'].to_numpy(), 0.0)
    np.testing.assert_array_equal(ctrl['auc'].to_numpy(), 0.5)

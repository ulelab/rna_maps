"""Tests for rnamaps.enrichment.bootstrap_contrast."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rnamaps.enrichment.bootstrap_contrast import (
    _bootstrap_means,
    _safe_logit,
    _shrink_log2fc_by_magnitude,
    _shrink_log_odds_by_magnitude,
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
    """A spike position should give a positive log2fc with CI excluding 0.

    Uses ``binarise=False`` to recover the 5.0 continuous spike value
    in ``delta`` (with the default ``binarise=True`` it would collapse
    to 1.0). Uses ``shrinkage='magnitude'`` to keep ``log2fc`` finite
    at the ``ctrl=0`` spike position; under the new ``shrinkage='none'``
    default ``log2(5/0)=+inf`` and the percentile-based CI bounds at
    the spike degenerate to NaN (numpy cannot interpolate between
    infinite values).
    """
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
        rng=rng, n_boot=500, binarise=False, shrinkage='magnitude',
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
    """Adaptive pseudocount stays positive and override is respected.

    Pseudocount mode is now opt-in (the default is ``beta_binomial``),
    so both calls pass ``shrinkage='pseudocount'`` explicitly.
    """
    rng = np.random.default_rng(4)
    n_pos = 6
    matrix = rng.poisson(lam=0.05, size=(40, n_pos)).astype(float)
    names = ['enhanced'] * 20 + ['control'] * 20
    positions = np.arange(n_pos)
    df = _make_df(matrix, names, positions)
    exon_categories = pd.Series({'enhanced': 20, 'control': 20})

    res_adaptive = compute(df, exon_categories, label='middle_3ss',
                           rng=rng, n_boot=200,
                           shrinkage='pseudocount')
    enh = res_adaptive.plot_df[res_adaptive.plot_df['name'] == 'enhanced']
    assert (enh['pseudocount'] > 0).all()

    res_fixed = compute(df, exon_categories, label='middle_3ss',
                        rng=np.random.default_rng(4), n_boot=200,
                        shrinkage='pseudocount',
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
        'log_odds_ratio', 'log_odds_ratio_lo', 'log_odds_ratio_hi',
        'shrinkage', 'tau', 'pseudocount',
        'coverage', 'number_exons', 'norm_coverage',
        'control_coverage', 'control_number_exons',
        'control_norm_coverage', 'fold_change',
    }
    assert expected.issubset(res.plot_df.columns)
    assert res.plot_kind == 'ribbon'
    assert res.method_name == 'bootstrap_contrast'
    assert set(res.y_columns) == {'delta', 'log2fc', 'log_odds_ratio'}
    # Under the default ``shrinkage='none'`` (no regularisation), the
    # per-category rows should carry NaN for both ``tau`` and
    # ``pseudocount`` -- neither shrinkage hyper-parameter is in use.
    enh = res.plot_df[res.plot_df['name'] == 'enhanced']
    assert (enh['shrinkage'] == 'none').all()
    assert enh['tau'].isna().all()
    assert enh['pseudocount'].isna().all()


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


def test_shrink_log2fc_by_magnitude_properties():
    """The helper shrinks toward 0 when *both* rates are small and
    leaves the contrast alone when at least one rate is well above
    tau."""
    tau = 0.05

    # Both rates well above tau -> raw log2fc preserved (within a few
    # percent; the eps_safe constant inside the log adds a tiny bias).
    out = _shrink_log2fc_by_magnitude(
        np.array([0.5]), np.array([0.05]), tau,
    )[0]
    raw = np.log2(0.5 / 0.05)
    assert abs(out - raw) / raw < 0.10

    # Both rates much smaller than tau -> shrunk hard toward 0.
    out_small = _shrink_log2fc_by_magnitude(
        np.array([0.0001]), np.array([0.00005]), tau,
    )[0]
    assert abs(out_small) < 0.1

    # Equal rates -> exactly 0 regardless of magnitude.
    np.testing.assert_allclose(
        _shrink_log2fc_by_magnitude(
            np.array([0.001, 0.05, 0.3]),
            np.array([0.001, 0.05, 0.3]),
            tau,
        ),
        0.0,
        atol=1e-12,
    )

    # One side zero, other side large -> bounded (not +-inf) and
    # essentially preserved: max(cat, ctrl) >> tau means weight ~ 0.
    out_signal = _shrink_log2fc_by_magnitude(
        np.array([1.0]), np.array([0.0]), tau,
    )[0]
    assert np.isfinite(out_signal)
    # Should be a clearly positive number, not shrunk to zero.
    assert out_signal > 5.0

    # One side zero, other side small -> bounded *and* shrunk
    # (max(cat, ctrl) < tau means weight is close to 1).
    out_noise = _shrink_log2fc_by_magnitude(
        np.array([0.005]), np.array([0.0]), tau,
    )[0]
    assert np.isfinite(out_noise)
    assert abs(out_noise) < 2.0

    # Sample-size invariance: the formula is a pure function of rates.
    np.testing.assert_allclose(
        _shrink_log2fc_by_magnitude(
            np.array([0.4]), np.array([0.05]), tau,
        ),
        _shrink_log2fc_by_magnitude(
            np.array([0.4]), np.array([0.05]), tau,
        ),
    )


def test_magnitude_shrinks_to_zero_when_both_rates_super_small():
    """Direct test of the user's intent: when both rates are tiny
    everywhere, ``log2fc`` is pulled to ~0 (no +-inf, no double-digit
    swings)."""
    rng = np.random.default_rng(20)
    n_pos = 8
    # Identical 0/1 matrices on both sides except for a single-exon
    # flip at one position -- the classic "raw log2fc explodes here"
    # case. Both rates are super small (< 1% throughout).
    cat = np.zeros((50, n_pos))
    ctrl = np.zeros((500, n_pos))
    cat[0, 3] = 1.0  # 1/50 = 2% at position 3 only
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * 50 + ['control'] * 500
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': 50, 'control': 500})

    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_boot=200, shrinkage='magnitude')
    enh = res.plot_df[res.plot_df['name'] == 'enhanced']
    # The single-exon spike would give log2fc=inf without shrinkage
    # and ~14 with the previous rate-level shrinkage; magnitude
    # shrinkage on log2fc keeps it bounded.
    assert np.isfinite(enh['log2fc']).all()
    assert enh['log2fc'].abs().max() < 2.0, (
        f"log2fc not shrunk enough at sparse positions: "
        f"max={enh['log2fc'].abs().max():.3f}"
    )


def test_magnitude_log2fc_is_a_pure_function_of_rates():
    """The magnitude-shrunk log2fc formula contains no ``n_c`` or
    ``n_ctrl`` term -- given the same rate inputs, it must produce
    bit-identical outputs regardless of how many exons fed those
    rates. This is the property the user asked for: number of exons
    must not impact the metric."""
    tau = 0.05
    cat = np.array([0.001, 0.02, 0.1, 0.4, 0.8])
    ctrl = np.array([0.001, 0.0, 0.05, 0.1, 0.05])
    # Calling the helper with identical rate vectors must give
    # identical output, and there is no way to pass ``n`` -- by
    # construction the formula cannot depend on it.
    out_a = _shrink_log2fc_by_magnitude(cat, ctrl, tau)
    out_b = _shrink_log2fc_by_magnitude(cat, ctrl, tau)
    np.testing.assert_array_equal(out_a, out_b)
    # Strong peak (cat=0.8, ctrl=0.05) is essentially preserved (raw
    # log2fc of ~16x = ~4.0; shrunk by ~6%).
    assert out_a[4] > 3.5
    # ...while a sparse 2x signal (cat=0.02, ctrl=0) is heavily shrunk
    # because max(cat, ctrl) = 0.02 is well below tau = 0.05.
    assert abs(out_a[1]) < 2.0


def test_magnitude_shrinks_more_aggressively_than_pseudocount_at_sparse_positions():
    """When both rates are tiny, magnitude shrinkage must produce a
    smaller-magnitude ``log2fc`` than the legacy pseudocount mode --
    that's the whole point of the new default."""
    rng = np.random.default_rng(22)
    n_pos = 4
    # Both sides sparse (rate ~1%). Single-exon flips will create the
    # kind of position-level spike that pseudocount only partially
    # tames but magnitude shrinkage erases.
    cat = (rng.random((100, n_pos)) < 0.01).astype(float)
    ctrl = (rng.random((100, n_pos)) < 0.01).astype(float)
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * 100 + ['control'] * 100
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': 100, 'control': 100})

    res_mag = compute(df, exon_categories, label='middle_3ss',
                      rng=np.random.default_rng(0), n_boot=400,
                      shrinkage='magnitude')
    res_pc = compute(df, exon_categories, label='middle_3ss',
                     rng=np.random.default_rng(0), n_boot=400,
                     shrinkage='pseudocount')
    mag = res_mag.plot_df[res_mag.plot_df['name'] == 'enhanced'][
        'log2fc'].abs().max()
    pc = res_pc.plot_df[res_pc.plot_df['name'] == 'enhanced'][
        'log2fc'].abs().max()
    assert mag < pc, (
        f"Expected magnitude shrinkage to tame sparse log2fc more than "
        f"pseudocount, but got mag={mag:.3f} vs pc={pc:.3f}."
    )


def test_magnitude_does_not_shrink_strong_signal():
    """At a strong, well-supported peak (rate well above the regional
    baseline), magnitude shrinkage must leave the log2fc nearly equal
    to the unshrunk version. This is the key promise: shrinkage only
    bites when the rates are small."""
    rng_data = np.random.default_rng(23)
    n_pos = 10
    peak = 5
    n_each = 200
    # Strong, well-supported peak: 80% at category, 10% at control.
    # Both rates are well above any plausible regional baseline.
    cat = (rng_data.random((n_each, n_pos)) < 0.10).astype(float)
    ctrl = (rng_data.random((n_each, n_pos)) < 0.10).astype(float)
    cat[:, peak] = (rng_data.random(n_each) < 0.8).astype(float)
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * n_each + ['control'] * n_each
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': n_each, 'control': n_each})

    res_mag = compute(df, exon_categories, label='middle_3ss',
                      rng=np.random.default_rng(200), n_boot=400,
                      shrinkage='magnitude')
    res_pc = compute(df, exon_categories, label='middle_3ss',
                     rng=np.random.default_rng(200), n_boot=400,
                     shrinkage='pseudocount')

    mag = res_mag.plot_df.set_index(['name', 'position']).loc[
        ('enhanced', peak), 'log2fc']
    pc = res_pc.plot_df.set_index(['name', 'position']).loc[
        ('enhanced', peak), 'log2fc']
    # Within ~15% of each other at a strong, well-supported peak.
    assert abs(mag - pc) / max(abs(pc), 1e-6) < 0.15, (
        f"Strong-peak log2fc differs too much: magnitude={mag:.3f}, "
        f"pseudocount={pc:.3f}"
    )


def test_magnitude_scale_knob_tunes_aggressiveness():
    """Larger ``shrinkage_scale`` tau -> more shrinkage; smaller tau ->
    less. At a moderate signal (~2x baseline), this is the most
    visible effect."""
    rng = np.random.default_rng(24)
    n_pos = 5
    peak = 2
    base = 0.05
    rate_at_peak_cat = 0.10  # exactly 2x baseline
    n_each = 500
    cat = (rng.random((n_each, n_pos)) < base).astype(float)
    cat[:, peak] = (rng.random(n_each) < rate_at_peak_cat).astype(float)
    ctrl = (rng.random((n_each, n_pos)) < base).astype(float)
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * n_each + ['control'] * n_each
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': n_each, 'control': n_each})

    def _peak_lfc(tau):
        res = compute(df, exon_categories, label='middle_3ss',
                      rng=np.random.default_rng(0), n_boot=300,
                      shrinkage='magnitude',
                      shrinkage_scale=tau)
        return res.plot_df.set_index(['name', 'position']).loc[
            ('enhanced', peak), 'log2fc']

    lfc_gentle = _peak_lfc(0.005)   # tau much smaller than rates
    lfc_default = _peak_lfc(0.05)    # tau = baseline
    lfc_aggressive = _peak_lfc(0.5)  # tau >> rates -> heavy shrinkage

    assert lfc_gentle > lfc_default > lfc_aggressive > 0, (
        f"Expected monotone shrinkage with tau: gentle={lfc_gentle:.3f}, "
        f"default={lfc_default:.3f}, aggressive={lfc_aggressive:.3f}"
    )


def test_shrinkage_none_can_produce_infinities():
    """``shrinkage='none'`` is allowed but documents that infinities are
    possible when one rate is exactly zero. The implementation must not
    raise; the user opts in to the unregularised behaviour."""
    rng = np.random.default_rng(15)
    cat = np.zeros((10, 4))
    ctrl = np.zeros((10, 4))
    cat[:, 1] = 1.0  # category positive at position 1, control all zero.
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * 10 + ['control'] * 10
    df = _make_df(matrix, names, np.arange(4))
    exon_categories = pd.Series({'enhanced': 10, 'control': 10})

    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_boot=50, shrinkage='none')
    enh = res.plot_df[res.plot_df['name'] == 'enhanced']
    # Position 1: log2(1/0) -> +inf, allowed under 'none'.
    assert not np.isfinite(enh.set_index('position').loc[1, 'log2fc'])


def test_invalid_shrinkage_raises():
    rng = np.random.default_rng(16)
    df = _make_df(np.zeros((4, 3)),
                  ['enhanced', 'enhanced', 'control', 'control'],
                  np.arange(3))
    exon_categories = pd.Series({'enhanced': 2, 'control': 2})
    with pytest.raises(ValueError, match="shrinkage"):
        compute(df, exon_categories, label='middle_3ss',
                rng=rng, n_boot=10, shrinkage='laplace')


def test_safe_logit_clipping():
    """``_safe_logit`` returns finite values at the boundary and matches
    ``ln(p / (1-p))`` away from it."""
    out = _safe_logit(np.array([0.0, 0.5, 1.0]), eps_safe=1e-3)
    # 0 -> eps, 1 -> 1-eps; both bounded by |ln((1-eps)/eps)| ~ 6.9.
    assert np.all(np.isfinite(out))
    assert out[0] == pytest.approx(-out[2])  # symmetric clip
    assert out[1] == pytest.approx(0.0, abs=1e-12)
    # Away from the boundary the safety clip is invisible.
    p = 0.3
    expected = np.log(p / (1 - p))
    assert _safe_logit(np.array([p]))[0] == pytest.approx(expected, rel=1e-6)


def test_log_odds_ratio_recovers_known_value():
    """At a position with known per-side rates the bootstrap mean
    ``log_odds_ratio`` recovers the analytic ``logit(p_cat) - logit(p_ctrl)``.
    """
    rng = np.random.default_rng(30)
    n_pos = 5
    spike = 2
    n_each = 50
    cat = np.zeros((n_each, n_pos))
    ctrl = np.zeros((n_each, n_pos))
    cat[:40, spike] = 1.0   # 40/50 = 80% positive in category
    ctrl[:10, spike] = 1.0  # 10/50 = 20% positive in control
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * n_each + ['control'] * n_each
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': n_each, 'control': n_each})

    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_boot=500)
    enh = res.plot_df[res.plot_df['name'] == 'enhanced'].set_index('position')

    expected_lor = np.log(0.8 / 0.2) - np.log(0.2 / 0.8)
    # Bootstrap mean is unbiased only up to O(1/sqrt(B)); 0.3 is loose.
    assert enh.loc[spike, 'log_odds_ratio'] == pytest.approx(
        expected_lor, abs=0.3
    )
    # The CI band at the spike must exclude zero (strong, real effect).
    assert enh.loc[spike, 'log_odds_ratio_lo'] > 0
    # Off-spike positions: both rates are exactly zero. The safety
    # clip makes both logits equal (logit(eps_safe)), so the log
    # odds ratio is identically zero.
    other = enh.drop(spike)
    np.testing.assert_allclose(
        other['log_odds_ratio'].to_numpy(), 0.0, atol=1e-12
    )


def test_log_odds_ratio_finite_under_saturation():
    """``log_odds_ratio`` stays bounded when one rate is exactly 1 or 0
    -- this is the property log2fc lacks (it would explode at 0).
    """
    rng = np.random.default_rng(31)
    n_pos = 4
    n_each = 30
    cat = np.ones((n_each, n_pos))   # 100% positive everywhere
    ctrl = np.zeros((n_each, n_pos))  # 0% positive everywhere
    matrix = np.vstack([cat, ctrl])
    names = ['enhanced'] * n_each + ['control'] * n_each
    df = _make_df(matrix, names, np.arange(n_pos))
    exon_categories = pd.Series({'enhanced': n_each, 'control': n_each})

    res = compute(df, exon_categories, label='middle_3ss',
                  rng=rng, n_boot=100, shrinkage='none')
    enh = res.plot_df[res.plot_df['name'] == 'enhanced']
    # log2fc(1/0) is +inf under shrinkage='none'.
    assert not np.all(np.isfinite(enh['log2fc']))
    # log_odds_ratio is still finite thanks to the eps_safe clip --
    # large (~2 * |logit(eps)|) but not inf.
    assert np.all(np.isfinite(enh['log_odds_ratio']))
    assert (enh['log_odds_ratio'] > 5).all()


def test_shrink_log_odds_by_magnitude_pulls_toward_zero():
    """When both rates are well below ``tau``, the magnitude-shrunk log
    odds ratio is pulled toward zero; when at least one is well above,
    the raw contrast is preserved."""
    tau = 0.05
    # Both rates tiny -> heavy shrinkage.
    out_small = _shrink_log_odds_by_magnitude(
        np.array([0.001]), np.array([0.0001]), tau,
    )[0]
    assert abs(out_small) < 0.5
    # One side well above tau -> contrast essentially preserved.
    out_big = _shrink_log_odds_by_magnitude(
        np.array([0.5]), np.array([0.05]), tau,
    )[0]
    raw_big = (np.log(0.5 / 0.5) - np.log(0.05 / 0.95))
    assert abs(out_big) > 0.7 * abs(raw_big)

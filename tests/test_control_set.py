"""Tests for rnamaps.preprocessing.apply_control_set."""

from __future__ import annotations

import pandas as pd
import pytest

from rnamaps.preprocessing import apply_control_set


def _make_df():
    rows = []
    # 5 strict-eligible controls.
    for i in range(5):
        rows.append({'category': 'control', 'dPSI': 0.005, 'FDR': 0.8})
    # 5 contaminated controls (large dPSI, low FDR).
    for i in range(5):
        rows.append({'category': 'control', 'dPSI': 0.04, 'FDR': 0.2})
    # 3 constitutive.
    for i in range(3):
        rows.append({'category': 'constitutive', 'dPSI': 0.01, 'FDR': 0.9})
    # 2 enhanced.
    for i in range(2):
        rows.append({'category': 'enhanced', 'dPSI': -0.2, 'FDR': 0.05})
    return pd.DataFrame(rows)


def test_default_is_noop():
    df = _make_df()
    out = apply_control_set(df, mode='default')
    pd.testing.assert_frame_equal(df.reset_index(drop=True),
                                  out.reset_index(drop=True))


def test_strict_filters_contaminated_controls():
    df = _make_df()
    out = apply_control_set(df, mode='strict',
                            control_max_dpsi=0.01,
                            control_min_fdr=0.5)
    n_ctrl = (out['category'] == 'control').sum()
    assert n_ctrl == 5
    # Other categories unchanged.
    assert (out['category'] == 'constitutive').sum() == 3
    assert (out['category'] == 'enhanced').sum() == 2


def test_constitutive_only_relabels_and_drops():
    df = _make_df()
    out = apply_control_set(df, mode='constitutive_only')
    assert (out['category'] == 'control').sum() == 3
    assert (out['category'] == 'constitutive').sum() == 0
    # Enhanced untouched.
    assert (out['category'] == 'enhanced').sum() == 2


def test_constitutive_only_raises_when_missing():
    df = _make_df()
    df = df[df['category'] != 'constitutive']
    with pytest.raises(ValueError, match="constitutive"):
        apply_control_set(df, mode='constitutive_only')


def test_unknown_mode_raises():
    df = _make_df()
    with pytest.raises(ValueError, match="Unknown"):
        apply_control_set(df, mode='not_a_mode')


def test_strict_skips_when_dpsi_missing():
    """VastDB-like input lacks dPSI; strict logs a warning and falls
    back to FDR-only filtering."""
    df = _make_df().drop(columns=['dPSI'])
    out = apply_control_set(df, mode='strict',
                            control_max_dpsi=0.01,
                            control_min_fdr=0.5)
    # FDR > 0.5 keeps the 5 strict-eligible controls and the 3 constitutive.
    n_ctrl = (out['category'] == 'control').sum()
    assert n_ctrl == 5

"""Tests for TPM attachment and expression-based control matching."""

from __future__ import annotations

import numpy as np
import pandas as pd

from rnamaps.expression_matching import (
    attach_tpm_to_exons,
    load_tpm_table,
    match_controls_by_expression,
)


def test_load_tpm_table_basic(tmp_path):
    tpm_path = tmp_path / "gene_tpm.tsv"
    tpm_path.write_text(
        "gene_id\ttpm\n"
        "ENSG000001\t10\n"
        "ENSG000001\t30\n"
        "TP53\t5\n"
        "BAD\tNA\n"
    )
    out = load_tpm_table(str(tpm_path))
    assert set(out.columns) == {"gene_id", "tpm"}
    assert len(out) == 2
    assert float(out.loc[out["gene_id"] == "ENSG000001", "tpm"].iloc[0]) == 20.0


def test_attach_picks_ensembl_when_better():
    df = pd.DataFrame(
        {
            "category": ["enhanced", "silenced", "control", "constitutive"],
            "GeneID": ["ENSG000001.2", "ENSG000002.5", "ENSG000003.1", "ENSG000004.3"],
            "geneSymbol": ["A1", "B2", "C3", "D4"],
        }
    )
    tpm_df = pd.DataFrame(
        {
            "gene_id": ["ENSG000001", "ENSG000002", "ENSG000003", "ENSG000004", "A1"],
            "tpm": [10.0, 20.0, 30.0, 40.0, 1.0],
        }
    )
    out, info = attach_tpm_to_exons(df, tpm_df)
    assert info["id_type"] == "ensembl"
    assert out["tpm"].notna().sum() == 4
    assert out["gene_id"].tolist() == [
        "ensg000001",
        "ensg000002",
        "ensg000003",
        "ensg000004",
    ]


def test_attach_picks_symbol_when_better():
    df = pd.DataFrame(
        {
            "category": ["enhanced", "silenced", "control", "constitutive"],
            "GeneID": ["ENSG999001.2", "ENSG999002.5", "ENSG999003.1", "ENSG999004.3"],
            "GENE": ["A1", "B2", "C3", "D4"],
        }
    )
    tpm_df = pd.DataFrame(
        {
            "gene_id": ["A1", "B2", "C3", "D4", "ENSG000001"],
            "tpm": [10.0, 20.0, 30.0, 40.0, 2.0],
        }
    )
    out, info = attach_tpm_to_exons(df, tpm_df)
    assert info["id_type"] == "symbol"
    assert out["tpm"].notna().sum() == 4
    assert out["gene_id"].tolist() == ["a1", "b2", "c3", "d4"]


def _make_expression_df():
    rows = []
    idx = 0
    # Regulated set: 5 low, 10 mid, 5 high TPM.
    for _ in range(5):
        rows.append({"idx": idx, "category": "enhanced", "tpm": 2.0})
        idx += 1
    for _ in range(10):
        rows.append({"idx": idx, "category": "silenced", "tpm": 10.0})
        idx += 1
    for _ in range(5):
        rows.append({"idx": idx, "category": "enhanced", "tpm": 60.0})
        idx += 1

    # Controls and constitutives include many extra points outside range.
    for _ in range(20):
        rows.append({"idx": idx, "category": "control", "tpm": 2.0})
        idx += 1
    for _ in range(20):
        rows.append({"idx": idx, "category": "control", "tpm": 10.0})
        idx += 1
    for _ in range(20):
        rows.append({"idx": idx, "category": "control", "tpm": 60.0})
        idx += 1
    for _ in range(5):
        rows.append({"idx": idx, "category": "control", "tpm": 300.0})
        idx += 1

    for _ in range(12):
        rows.append({"idx": idx, "category": "constitutive", "tpm": 2.0})
        idx += 1
    for _ in range(12):
        rows.append({"idx": idx, "category": "constitutive", "tpm": 10.0})
        idx += 1
    for _ in range(12):
        rows.append({"idx": idx, "category": "constitutive", "tpm": 60.0})
        idx += 1
    for _ in range(4):
        rows.append({"idx": idx, "category": "constitutive", "tpm": 300.0})
        idx += 1

    df = pd.DataFrame(rows).set_index("idx")
    df["gene_id"] = [f"gene{i}" for i in range(len(df))]
    return df


def _bin_counts(values, edges):
    bins = pd.cut(values, bins=edges, labels=False, include_lowest=True, right=True)
    return bins.value_counts().sort_index()


def test_match_preserves_distribution():
    df = _make_expression_df()
    out, summary = match_controls_by_expression(
        df,
        n_bins=3,
        pseudocount=1.0,
        min_tpm=0.0,
        also_constitutive=True,
        rng=np.random.default_rng(42),
    )
    edges = summary["regulated_bin_edges"]
    reg = out[out["category"].isin({"enhanced", "silenced"})]
    ctrl = out[out["category"] == "control"]
    reg_counts = _bin_counts(np.log10(reg["tpm"] + 1.0), edges)
    ctrl_counts = _bin_counts(np.log10(ctrl["tpm"] + 1.0), edges)
    reg_fracs = reg_counts / reg_counts.sum()
    expected_ctrl = reg_fracs * ctrl_counts.sum()
    assert all(abs(ctrl_counts - expected_ctrl).fillna(0) <= 1.0)


def test_match_drops_constitutive_when_optout():
    df = _make_expression_df()
    before_const = int((df["category"] == "constitutive").sum())
    out, _ = match_controls_by_expression(
        df,
        n_bins=3,
        pseudocount=1.0,
        min_tpm=0.0,
        also_constitutive=False,
        rng=np.random.default_rng(1),
    )
    after_const = int((out["category"] == "constitutive").sum())
    assert before_const == after_const


def test_match_handles_empty_negative_bin():
    df = pd.DataFrame(
        {
            "category": ["enhanced"] * 5 + ["silenced"] * 5 + ["control"] * 6,
            "tpm": [2.0] * 5 + [50.0] * 5 + [2.0] * 6,
            "gene_id": [f"g{i}" for i in range(16)],
        }
    )
    out, _ = match_controls_by_expression(
        df,
        n_bins=2,
        pseudocount=1.0,
        min_tpm=0.0,
        also_constitutive=False,
        rng=np.random.default_rng(0),
    )
    assert (out["category"] == "control").sum() == 0


def test_match_seeded_reproducible():
    df = _make_expression_df()
    out1, _ = match_controls_by_expression(
        df,
        n_bins=3,
        pseudocount=1.0,
        min_tpm=0.0,
        also_constitutive=True,
        rng=np.random.default_rng(123),
    )
    out2, _ = match_controls_by_expression(
        df,
        n_bins=3,
        pseudocount=1.0,
        min_tpm=0.0,
        also_constitutive=True,
        rng=np.random.default_rng(123),
    )
    assert sorted(out1.index.tolist()) == sorted(out2.index.tolist())

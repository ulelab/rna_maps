"""Gene-expression aware matching of negative exon sets."""

from __future__ import annotations

import logging
import re
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd


REGULATED_CATEGORIES = {"enhanced", "silenced"}
NEGATIVE_CATEGORIES = ("control", "constitutive")


def _normalise_ensembl(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().strip('"').strip("'")
    if not text:
        return None
    return re.sub(r"\.\d+$", "", text).lower()


def _normalise_symbol(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().strip('"').strip("'")
    if not text:
        return None
    return text.lower()


def _series_from_first_present(df: pd.DataFrame, columns: Iterable[str]) -> pd.Series:
    for col in columns:
        if col in df.columns:
            return df[col]
    return pd.Series(index=df.index, dtype=object)


def load_tpm_table(path: str) -> pd.DataFrame:
    """Load a gene TPM table with columns ``gene_id`` and ``tpm``."""
    df = pd.read_csv(path, sep=None, engine="python", comment="#")
    if df.empty:
        raise ValueError(f"TPM table is empty: {path}")

    df.columns = [str(c).strip().lower() for c in df.columns]
    if "gene_id" in df.columns and "tpm" in df.columns:
        parsed = df.loc[:, ["gene_id", "tpm"]].copy()
    elif len(df.columns) >= 2:
        parsed = df.iloc[:, :2].copy()
        parsed.columns = ["gene_id", "tpm"]
        logging.warning(
            "TPM table lacks explicit 'gene_id'/'tpm' headers; using first "
            "two columns."
        )
    else:
        raise ValueError(
            "TPM table must have at least two columns: gene_id and tpm."
        )

    parsed["gene_id"] = parsed["gene_id"].astype(str).str.strip().str.strip('"').str.strip("'")
    parsed["tpm"] = pd.to_numeric(parsed["tpm"], errors="coerce")
    parsed = parsed.dropna(subset=["gene_id", "tpm"])
    parsed = parsed[parsed["gene_id"] != ""]
    parsed = parsed[parsed["tpm"] >= 0]

    if parsed.empty:
        raise ValueError(f"No valid gene_id/tpm rows found in: {path}")

    # Duplicated gene identifiers are averaged as a pragmatic fallback.
    parsed = parsed.groupby("gene_id", as_index=False, sort=False)["tpm"].mean()
    logging.info(f"Loaded TPM table with {len(parsed)} unique gene IDs from {path}")
    return parsed


def attach_tpm_to_exons(
    df_exons: pd.DataFrame, tpm_df: pd.DataFrame, min_match_rate: float = 0.5
) -> Tuple[pd.DataFrame, Dict[str, float | str | int]]:
    """Attach TPM to exon rows, auto-detecting Ensembl vs symbol IDs."""
    df = df_exons.copy()

    ens_series = _series_from_first_present(df, ["GeneID", "gene_id"]).map(_normalise_ensembl)
    sym_series = _series_from_first_present(df, ["GENE", "geneSymbol", "gene_name"]).map(
        _normalise_symbol
    )

    tpm_ens = tpm_df.copy()
    tpm_ens["key"] = tpm_ens["gene_id"].map(_normalise_ensembl)
    tpm_ens = tpm_ens.dropna(subset=["key"]).groupby("key", as_index=False)["tpm"].mean()

    tpm_sym = tpm_df.copy()
    tpm_sym["key"] = tpm_sym["gene_id"].map(_normalise_symbol)
    tpm_sym = tpm_sym.dropna(subset=["key"]).groupby("key", as_index=False)["tpm"].mean()

    map_ens = dict(zip(tpm_ens["key"], tpm_ens["tpm"]))
    map_sym = dict(zip(tpm_sym["key"], tpm_sym["tpm"]))

    ens_valid = ens_series.notna().sum()
    sym_valid = sym_series.notna().sum()
    ens_hits = ens_series.map(map_ens).notna().sum()
    sym_hits = sym_series.map(map_sym).notna().sum()
    ens_rate = (ens_hits / ens_valid) if ens_valid else 0.0
    sym_rate = (sym_hits / sym_valid) if sym_valid else 0.0

    if ens_rate >= sym_rate:
        chosen = "ensembl"
        df["gene_id"] = ens_series
        df["tpm"] = ens_series.map(map_ens)
        chosen_rate = ens_rate
    else:
        chosen = "symbol"
        df["gene_id"] = sym_series
        df["tpm"] = sym_series.map(map_sym)
        chosen_rate = sym_rate

    logging.info(
        "TPM ID auto-detect: ensembl match %.1f%% (%s/%s), symbol match %.1f%% (%s/%s); "
        "using %s IDs.",
        ens_rate * 100.0,
        ens_hits,
        ens_valid,
        sym_rate * 100.0,
        sym_hits,
        sym_valid,
        chosen,
    )

    if chosen_rate < min_match_rate:
        raise ValueError(
            "TPM ID auto-detect failed: low overlap between exons and TPM table "
            f"(ensembl={ens_rate:.1%}, symbol={sym_rate:.1%})."
        )

    missing_by_cat = (
        df.assign(_missing_tpm=df["tpm"].isna())
        .groupby("category", dropna=False)["_missing_tpm"]
        .sum()
        .to_dict()
    )
    if any(v > 0 for v in missing_by_cat.values()):
        logging.warning("Rows with missing TPM after ID matching by category: %s", missing_by_cat)

    info: Dict[str, float | str | int] = {
        "id_type": chosen,
        "ensembl_rate": float(ens_rate),
        "symbol_rate": float(sym_rate),
        "attached_rows": int(df["tpm"].notna().sum()),
        "total_rows": int(len(df)),
    }
    return df, info


def _allocate_counts(fracs: np.ndarray, total: int, avail: np.ndarray) -> np.ndarray:
    desired = fracs * float(total)
    counts = np.floor(desired).astype(int)
    counts = np.minimum(counts, avail)
    remaining = int(total - counts.sum())
    if remaining <= 0:
        return counts

    # Largest remainder allocation while respecting per-bin availability.
    order = np.argsort(-(desired - np.floor(desired)))
    for idx in order:
        if remaining == 0:
            break
        spare = int(avail[idx] - counts[idx])
        if spare <= 0:
            continue
        add = min(spare, remaining)
        counts[idx] += add
        remaining -= add
    return counts


def _subsample_negative_category(
    df_neg: pd.DataFrame, reg_fracs: np.ndarray, n_bins: int, rng: np.random.Generator
) -> pd.DataFrame:
    avail = (
        df_neg["expr_bin"]
        .value_counts()
        .reindex(range(n_bins), fill_value=0)
        .to_numpy(dtype=int)
    )
    if reg_fracs.sum() == 0:
        return df_neg.iloc[0:0].copy()

    feasible = np.where(reg_fracs > 0, avail / reg_fracs, np.inf)
    n_keep = int(np.floor(np.min(feasible)))
    if n_keep <= 0:
        return df_neg.iloc[0:0].copy()

    target = _allocate_counts(reg_fracs, n_keep, avail)
    keep_idx = []
    for bin_idx in range(n_bins):
        n_take = int(target[bin_idx])
        if n_take <= 0:
            continue
        candidates = df_neg.index[df_neg["expr_bin"] == bin_idx].to_numpy()
        chosen = rng.choice(candidates, size=n_take, replace=False)
        keep_idx.extend(chosen.tolist())
    return df_neg.loc[keep_idx].copy()


def match_controls_by_expression(
    df_exons: pd.DataFrame,
    n_bins: int,
    pseudocount: float,
    min_tpm: float,
    also_constitutive: bool,
    rng: np.random.Generator,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Stratified expression matching for control (and optionally constitutive)."""
    if n_bins < 2:
        raise ValueError("--tpm_n_bins must be >= 2")
    if pseudocount <= 0:
        raise ValueError("--tpm_pseudocount must be > 0")

    df_all = df_exons.copy()
    before_counts = df_all["category"].value_counts().to_dict()

    exempt_frames = []
    if not also_constitutive:
        exempt_const = df_all[df_all["category"] == "constitutive"].copy()
        if not exempt_const.empty:
            exempt_frames.append(exempt_const)
        df = df_all[df_all["category"] != "constitutive"].copy()
    else:
        df = df_all

    valid_tpm_mask = df["tpm"].notna() & (df["tpm"] >= min_tpm)
    dropped = df.loc[~valid_tpm_mask, "category"].value_counts().to_dict()
    if dropped:
        logging.info(
            "Dropping rows with missing/low TPM (min_tpm=%s) by category: %s",
            min_tpm,
            dropped,
        )
    df = df.loc[valid_tpm_mask].copy()

    if df.empty:
        raise ValueError("No rows remain after TPM filtering.")

    df["log_tpm"] = np.log10(df["tpm"] + pseudocount)
    reg_mask = df["category"].isin(REGULATED_CATEGORIES)
    reg_df = df.loc[reg_mask]
    if reg_df.empty:
        raise ValueError("No regulated exons with valid TPM found for expression matching.")

    q = np.linspace(0, 1, n_bins + 1)
    edges = np.quantile(reg_df["log_tpm"].to_numpy(), q=q)
    edges = np.unique(edges)
    if len(edges) < 2:
        raise ValueError("Regulated-set TPM distribution has insufficient variation for binning.")
    n_bins_eff = len(edges) - 1

    df["expr_bin"] = pd.cut(
        df["log_tpm"], bins=edges, labels=False, include_lowest=True, right=True
    )
    df = df.dropna(subset=["expr_bin"]).copy()
    df["expr_bin"] = df["expr_bin"].astype(int)

    reg_df = df[df["category"].isin(REGULATED_CATEGORIES)]
    reg_counts = (
        reg_df["expr_bin"].value_counts().reindex(range(n_bins_eff), fill_value=0).to_numpy(dtype=int)
    )
    reg_total = int(reg_counts.sum())
    reg_fracs = reg_counts / reg_total

    negative_targets = ["control"]
    if also_constitutive:
        negative_targets.append("constitutive")

    keep_frames = [df[~df["category"].isin(negative_targets)].copy()]
    neg_summary = {}
    for cat in negative_targets:
        cat_df = df[df["category"] == cat].copy()
        if cat_df.empty:
            neg_summary[cat] = {"before": 0, "after": 0}
            continue
        kept = _subsample_negative_category(cat_df, reg_fracs, n_bins_eff, rng)
        keep_frames.append(kept)
        neg_summary[cat] = {"before": int(len(cat_df)), "after": int(len(kept))}
        if len(cat_df) > 0 and len(kept) / len(cat_df) < 0.2:
            logging.warning(
                "Expression matching removed >80%% of %s exons (%s -> %s).",
                cat,
                len(cat_df),
                len(kept),
            )

    out = pd.concat(keep_frames + exempt_frames, axis=0).sort_index()
    after_counts = out["category"].value_counts().to_dict()

    summary: Dict[str, object] = {
        "before_counts": {k: int(v) for k, v in before_counts.items()},
        "after_counts": {k: int(v) for k, v in after_counts.items()},
        "dropped_missing_or_low_tpm": {k: int(v) for k, v in dropped.items()},
        "regulated_bin_edges": [float(x) for x in edges],
        "regulated_bin_counts": [int(x) for x in reg_counts.tolist()],
        "negative_summary": neg_summary,
    }
    return out, summary

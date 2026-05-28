"""Per-bin true-positive / false-positive table for the cassette splice sites.

For each non-overlapping bin tiling the ``middle_3ss`` and ``middle_5ss``
windows, we compute the fraction of regulated exons with any CLIP signal
in the bin (TPR / sensitivity) and the analogous fraction in the control
set (FPR). Additional columns: enrichment (TPR/FPR), log2 enrichment,
two-sided Fisher's exact p-value on the 2x2 of ``{category, control} x
{bound, unbound}``, and a Benjamini-Hochberg q-value computed within
each ``(region, category)`` family.

This method writes a tabular output only -- it does not contribute to
the RNA map PDFs. Use it when you want a single interpretable table per
BED file describing where binding actually distinguishes regulated from
control exons at coarse (typically 50 nt) resolution.

Other splice-site regions (``upstream_5ss``, ``downstream_3ss``, etc.)
are skipped: by construction this analysis is about the cassette flanks.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats as _sstats

from rnamaps.enrichment._base import EnrichmentResult


_MIDDLE_LABELS = ("middle_3ss", "middle_5ss")


def _bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values for a vector of p-values.

    NaN p-values pass through as NaN and are excluded from the rank count.
    """
    pvals = np.asarray(pvals, dtype=np.float64)
    out = np.full_like(pvals, np.nan, dtype=np.float64)
    valid = ~np.isnan(pvals)
    if not valid.any():
        return out
    p = pvals[valid]
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, dtype=np.int64)
    ranks[order] = np.arange(1, n + 1)
    q = p * n / ranks
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_unsorted = np.empty(n, dtype=np.float64)
    q_unsorted[order] = q_sorted
    out[valid] = np.clip(q_unsorted, 0.0, 1.0)
    return out


def _bin_edges(n_positions: int, bin_size: int):
    """Splice-site-anchored 0-indexed column slices.

    Tiles upstream and downstream of the splice-site column
    independently so every bin lies entirely on one side of the SS and
    bin midpoints align to multiples of ``bin_size`` away from the SS.
    Tail bins narrower than ``bin_size`` are dropped (typically a single
    base when ``window`` isn't a clean multiple of ``bin_size``) so the
    output stays interpretable as full-width bins.

    Returns a list of ``(start, end)`` slices in column order
    (upstream-most → downstream-most).
    """
    # n_positions = 2*window + 1; SS sits at column index `window`.
    ss_idx = (n_positions - 1) // 2
    # Upstream: build outward from the SS so the bin immediately left of
    # the SS is full-width; then reverse to column order.
    up = []
    s = ss_idx
    while s - bin_size >= 0:
        up.append((s - bin_size, s))
        s -= bin_size
    up.reverse()
    # Downstream: the SS column itself goes into the first downstream
    # bin so the boundary aligns with the SS coordinate.
    down = []
    s = ss_idx
    while s + bin_size <= n_positions:
        down.append((s, s + bin_size))
        s += bin_size
    return up + down


def compute(
    region,
    exon_categories: pd.Series,
    label: str,
    *,
    bin_size: int = 50,
    control_label: str = "control",
) -> Optional[EnrichmentResult]:
    """Per-bin TPR / FPR / enrichment table for one splice-site region.

    Returns ``None`` for non-cassette regions so the caller can skip
    them transparently.

    Parameters
    ----------
    region : rnamaps.coverage.RegionCoverage
    exon_categories : pd.Series
    label : str
        Splice-site region label. Must be ``"middle_3ss"`` or
        ``"middle_5ss"`` for any work to happen.
    bin_size : int
        Bin width in nucleotides. The position vector is 1 nt resolution
        so this is also the column-stride for the bin.
    control_label : str
    """
    if label not in _MIDDLE_LABELS:
        return None
    if control_label not in exon_categories.index:
        return None

    matrix = region.matrix
    if matrix.size == 0:
        return None

    positions = region.positions
    n_pos = matrix.shape[1]
    bins = _bin_edges(n_pos, bin_size)
    if not bins:
        return None

    # Splice-site coordinate: positions are 1..n_pos with the SS at the
    # centre column after slop. ``window`` is (n_pos - 1) / 2.
    ss_position = (int(positions[0]) + int(positions[-1])) / 2.0

    # Collapse each bin to a per-exon 0/1 "any signal in this bin" mask
    # in one pass over the matrix.
    presence_per_bin = np.empty((matrix.shape[0], len(bins)), dtype=bool)
    for i, (s, e) in enumerate(bins):
        presence_per_bin[:, i] = (matrix[:, s:e] > 0).any(axis=1)

    ctrl_mask = (region.exon_names == control_label)
    n_ctrl = int(ctrl_mask.sum())
    if n_ctrl == 0:
        return None

    categories = [c for c in exon_categories.index if c != control_label]

    rows = []
    for cat in categories:
        cat_mask = (region.exon_names == cat)
        n_cat = int(cat_mask.sum())
        if n_cat == 0:
            continue
        n_cat_bound = presence_per_bin[cat_mask].sum(axis=0).astype(int)
        n_ctrl_bound = presence_per_bin[ctrl_mask].sum(axis=0).astype(int)
        for i, (s, e) in enumerate(bins):
            bin_start_pos = int(positions[s])
            bin_end_pos = int(positions[e - 1])
            bin_mid = (bin_start_pos + bin_end_pos) / 2.0
            n_cb = int(n_cat_bound[i])
            n_kb = int(n_ctrl_bound[i])
            tpr = n_cb / n_cat
            fpr = n_kb / n_ctrl
            if fpr > 0 and tpr > 0:
                enrichment = tpr / fpr
                log2_enrichment = float(np.log2(enrichment))
            elif fpr == 0 and tpr > 0:
                enrichment = np.inf
                log2_enrichment = np.inf
            elif tpr == 0 and fpr > 0:
                enrichment = 0.0
                log2_enrichment = -np.inf
            else:
                enrichment = np.nan
                log2_enrichment = np.nan
            try:
                _, pvalue = _sstats.fisher_exact(
                    [[n_cb, n_cat - n_cb],
                     [n_kb, n_ctrl - n_kb]],
                    alternative="two-sided",
                )
            except ValueError:
                pvalue = np.nan
            rows.append({
                "region": label,
                "category": cat,
                "bin_index": i,
                "bin_start": bin_start_pos,
                "bin_end": bin_end_pos,
                "bin_center_offset": float(bin_mid - ss_position),
                "bin_width": e - s,
                "n_regulated": n_cat,
                "n_regulated_bound": n_cb,
                "tpr": float(tpr),
                "n_control": n_ctrl,
                "n_control_bound": n_kb,
                "fpr": float(fpr),
                "enrichment": float(enrichment) if np.isfinite(enrichment)
                              else enrichment,
                "log2_enrichment": log2_enrichment,
                "pvalue": float(pvalue) if pvalue == pvalue else np.nan,
            })

    if not rows:
        return None

    df = pd.DataFrame(rows)
    df["qvalue"] = np.nan
    for (_, _), grp in df.groupby(["region", "category"], sort=False):
        df.loc[grp.index, "qvalue"] = _bh_qvalues(grp["pvalue"].to_numpy())

    return EnrichmentResult(
        plot_df=df,
        clusters_df=pd.DataFrame(),
        plot_kind="table_only",
        y_columns=(),
        method_name="bin_tpr",
        ylabel="",
    )

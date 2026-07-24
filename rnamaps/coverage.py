"""CLIP coverage calculation around splice sites.

The previous implementation invoked ``bedtools coverage -d`` (per-base
coverage) and rehydrated the full ``(n_exons × n_positions)`` long-form
DataFrame in pandas, then pivoted that back to a matrix downstream. For
realistic inputs (~80K exons × 601 positions × 4 regions ≈ 200M rows
with object-dtype columns) that approach peaked at 25-35 GB of RAM and
ran into multi-hour bedtools-d output parsing.

This module now:

- Uses ``bedtools intersect -wa -wb`` to get only the actual overlap
  pairs (one row per (exon-window, crosslink) hit) rather than emitting
  one row per (exon, base).
- Builds a compact ``(n_exons × n_positions)`` ``int32`` count matrix
  directly with vectorised ``numpy.add.at``.
- Exposes a small ``RegionCoverage`` value object so downstream
  enrichment / heatmap code can operate on the matrix directly without
  re-pivoting a 200M-row DataFrame.

The aggregated, category-level ``df_plot`` (used by the legacy Fisher
linegraph and the TSV outputs) is unchanged in shape and column names.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import pybedtools as pbt
import scipy.stats as stats


_VALID_SCORE_MODES = (
    "ignore", "raw", "per_transcript_zscore",
)


@dataclass
class RegionCoverage:
    """Per-exon × per-position coverage for one splice-site region.

    Attributes
    ----------
    label : str
        Splice-site region label, e.g. ``"middle_3ss"``.
    matrix : ndarray, shape (n_exons, n_positions)
        Raw integer crosslink counts. Downstream enrichment methods
        binarise this on the fly (``matrix > 0``) when they need the
        "fraction of exons positive" semantics; the legacy aggregated
        Fisher path uses the raw sums via :meth:`category_sums`.
    exon_ids : ndarray of str, shape (n_exons,)
        Exon identifier (without the leading ``category_`` prefix).
    exon_names : ndarray of str, shape (n_exons,)
        Per-exon category label (``enhanced``, ``silenced``, ``control``,
        ``constitutive``, ...).
    positions : ndarray of int, shape (n_positions,)
        1-indexed transcript-coordinate positions (1..2*window+1).
    """

    label: str
    matrix: np.ndarray
    exon_ids: np.ndarray
    exon_names: np.ndarray
    positions: np.ndarray

    @property
    def n_exons(self) -> int:
        return int(self.matrix.shape[0])

    @property
    def n_positions(self) -> int:
        return int(self.matrix.shape[1])

    def category_sums(self) -> dict:
        """Sum crosslink counts per category, per position.

        Returns
        -------
        dict[str, ndarray]
            Maps category name → ``(n_positions,)`` ``int64`` sum vector.
        """
        out = {}
        for cat in np.unique(self.exon_names):
            out[cat] = self.matrix[self.exon_names == cat].sum(axis=0)
        return out

    def select(self, categories):
        """Return ``(matrix, names)`` restricted to rows in ``categories``.

        Parameters
        ----------
        categories : iterable of str

        Returns
        -------
        sub_matrix : ndarray
        sub_names : ndarray of str
        """
        mask = np.isin(self.exon_names, list(categories))
        return self.matrix[mask], self.exon_names[mask]

    def save_npz(self, path: str) -> None:
        """Persist this RegionCoverage to ``path`` as a compressed .npz.

        Layout: ``matrix`` (n_exons × n_positions, raw integer counts before
        binarisation/smoothing), ``exon_ids``, ``exon_names``, ``positions``,
        ``label`` (0-d string array). Round-trippable via :meth:`load_npz`.
        """
        np.savez_compressed(
            path,
            matrix=self.matrix,
            exon_ids=np.asarray(self.exon_ids, dtype=str),
            exon_names=np.asarray(self.exon_names, dtype=str),
            positions=self.positions,
            label=np.array(self.label),
        )

    @classmethod
    def load_npz(cls, path: str) -> "RegionCoverage":
        """Inverse of :meth:`save_npz`."""
        with np.load(path, allow_pickle=False) as data:
            return cls(
                label=str(data['label']),
                matrix=data['matrix'],
                exon_ids=data['exon_ids'],
                exon_names=data['exon_names'],
                positions=data['positions'],
            )


def _ss_df_to_pbt(ss_df: pd.DataFrame, window: int, fai: str):
    """Build a sorted, windowed pybedtools BedTool from a splice-site frame.

    The input ``ss_df`` is expected to have BED6 columns
    (``chr, start, end, name, score, strand``); ``name`` carries the
    ``category_exon-id`` string built by :func:`get_ss_bed`. Rows with
    missing or placeholder names are dropped.
    """
    ss_df = ss_df.loc[
        ss_df['name'].notna() & (ss_df['name'].astype(str) != ".")
    ].copy()
    ss_df['name'] = ss_df['name'].astype(str)
    bed = pbt.BedTool.from_dataframe(
        ss_df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    ).sort().slop(l=window, r=window, s=True, g=fai)
    return ss_df, bed


def _overlap_pairs_to_matrix(
    overlap_bed: pbt.BedTool,
    name_to_row: dict,
    n_exons: int,
    n_pos: int,
    score_mode: str = "ignore",
) -> np.ndarray:
    """Convert bedtools intersect -wa -wb output to a (signal) matrix.

    Each output line carries A's BED6 (windowed exon) followed by B's
    fields (crosslink / score). We extract:

    - col 1 (A.start, 0-based start of the windowed exon)
    - col 2 (A.end, 0-based-exclusive end)
    - col 3 (A.name, ``category_exon-id``)
    - col 5 (A.strand)
    - col 7 (B.start, 0-based start of the B feature)
    - col 8 (B.end, 0-based-exclusive end of the B feature)
    - col 9 (B.name, the BED-track entry's name; used as the transcript
      identifier for the per-transcript score modes only)
    - col 10 (B.score, the BED-track entry's score; used when
      ``score_mode != 'ignore'``)

    Each B-feature deposits one count (or one ``b_score``, in score modes)
    at *every* base it covers inside the A-window — B is clipped to A,
    then expanded to per-base entries. For single-nucleotide iCLIP
    crosslink BEDs (``b_end == b_start + 1``) this collapses to the
    legacy single-position behaviour; for wider B features (e.g.
    finemapped windows, peaks) the signal is spread across all bases
    they cover. In the ``raw`` / per-transcript score modes each base
    receives the *same* score (the per-feature value is replicated, not
    divided by length); if you want length-normalised semantics, divide
    the BED score by feature length before passing it in.

    Position within the window is converted to transcript coordinates
    (1-indexed, 5'→3' of the exon strand) so + and - strand exons share
    the same coordinate system in the matrix.

    Parameters
    ----------
    score_mode : {"ignore", "raw", "per_transcript_zscore"}
        ``"ignore"`` (default) counts each overlap as ``+1`` (matrix
        dtype int32; exact legacy behaviour). ``"raw"`` uses the BED
        score column as-is (matrix becomes float64).
        ``"per_transcript_zscore"`` accumulates ``raw`` first, then
        renormalises rows of crosslinks sharing the same B.name
        (transcript) before the usual exon-row accumulation, producing
        zero-mean, unit-variance scores per transcript.
    """
    if score_mode not in _VALID_SCORE_MODES:
        raise ValueError(
            f"score_mode must be one of {_VALID_SCORE_MODES}, got "
            f"{score_mode!r}"
        )

    use_scores = score_mode != "ignore"
    dtype = np.float64 if use_scores else np.int32
    matrix = np.zeros((n_exons, n_pos), dtype=dtype)

    if not overlap_bed.fn or not _bed_has_content(overlap_bed.fn):
        return matrix

    # Read only the columns we need. bedtools intersect -wa -wb prefixes
    # A's columns then appends B's; A=0..5, B=6.. -- B may be BED3+ but
    # for the score modes we additionally need B.score (col 10) and
    # B.name (col 9, per-transcript modes only). We always read B.end
    # (col 8) so that multi-base B features get spread across all the
    # window bases they cover, not just their leftmost base.
    base_usecols = [1, 2, 3, 5, 7, 8]
    base_names = ['a_start', 'a_end', 'a_name', 'a_strand',
                  'b_start', 'b_end']
    base_dtype = {'a_start': np.int64, 'a_end': np.int64,
                  'a_name': str, 'a_strand': str,
                  'b_start': np.int64, 'b_end': np.int64}
    if score_mode == "raw":
        usecols = base_usecols + [10]
        names = base_names + ['b_score']
        dtype_map = {**base_dtype, 'b_score': np.float64}
    elif score_mode == "per_transcript_zscore":
        usecols = base_usecols + [9, 10]
        names = base_names + ['b_name', 'b_score']
        dtype_map = {**base_dtype, 'b_name': str, 'b_score': np.float64}
    else:
        usecols = base_usecols
        names = base_names
        dtype_map = base_dtype

    df = pd.read_csv(
        overlap_bed.fn,
        sep='\t',
        header=None,
        usecols=usecols,
        names=names,
        dtype=dtype_map,
    )
    if df.empty:
        return matrix

    if score_mode == "per_transcript_zscore":
        # Renormalise B.score *within each B.name (transcript)* across
        # ALL the entries that survived the intersect for this region.
        # Rows whose B.name is the BED placeholder "." are treated as
        # their own singleton transcripts (no normalisation possible).
        b_name = df['b_name'].astype(str).to_numpy()
        b_score = df['b_score'].to_numpy()
        renormed = b_score.copy()
        # Group transcripts (mask out singletons).
        is_singleton = (b_name == ".") | (b_name == "")
        if (~is_singleton).any():
            g = pd.Series(b_score[~is_singleton]).groupby(
                pd.Series(b_name[~is_singleton]).values
            )
            # ddof=0 to keep it well-defined at single-entry txs.
            renormed_vals = g.transform(
                lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else 0.0
            ).to_numpy()
            tmp = renormed.copy()
            tmp[~is_singleton] = renormed_vals
            renormed = tmp
        df = df.assign(b_score=renormed)

    # Map exon name -> row index (drop overlaps for names not present,
    # which shouldn't happen in practice but is cheap to guard against).
    rows = df['a_name'].map(name_to_row).to_numpy()
    valid_row = np.isfinite(rows.astype(np.float64, copy=False))
    if not valid_row.any():
        return matrix

    a_strand = df['a_strand'].to_numpy()[valid_row]
    a_start = df['a_start'].to_numpy()[valid_row]
    a_end = df['a_end'].to_numpy()[valid_row]
    b_start = df['b_start'].to_numpy()[valid_row]
    b_end = df['b_end'].to_numpy()[valid_row]
    rows = rows[valid_row].astype(np.int64, copy=False)
    if use_scores:
        scores = df['b_score'].to_numpy()[valid_row].astype(
            np.float64, copy=False
        )

    # Clip B to A: bedtools intersect -wa -wb reports the *original* B
    # interval, which may extend outside the A-window. We only want the
    # bases that actually fall inside the window.
    lo = np.maximum(b_start, a_start)
    hi = np.minimum(b_end, a_end)  # exclusive
    lengths = hi - lo
    nonempty = lengths > 0
    if not nonempty.any():
        return matrix

    lo = lo[nonempty]
    lengths = lengths[nonempty]
    a_strand = a_strand[nonempty]
    a_start = a_start[nonempty]
    a_end = a_end[nonempty]
    rows = rows[nonempty]
    if use_scores:
        scores = scores[nonempty]

    # Expand each overlap row into one entry per base it covers in the
    # window. For a 75-nt feature this produces 75 entries with the same
    # row index (and the same b_score, in score modes); for the 1-nt
    # iCLIP case this is a no-op repeat.
    rows_exp = np.repeat(rows, lengths)
    strand_exp = np.repeat(a_strand, lengths)
    a_start_exp = np.repeat(a_start, lengths)
    a_end_exp = np.repeat(a_end, lengths)
    # within-feature offset 0..length-1
    total = int(lengths.sum())
    group_starts = np.empty_like(lengths)
    group_starts[0] = 0
    if lengths.size > 1:
        group_starts[1:] = np.cumsum(lengths[:-1])
    within = np.arange(total, dtype=np.int64) - np.repeat(
        group_starts, lengths
    )
    genomic = np.repeat(lo, lengths) + within

    # Transcript-coordinate position (1-indexed):
    #   + strand: pos = genomic - a_start + 1
    #   - strand: pos = a_end - genomic
    pos = np.where(
        strand_exp == '+',
        genomic - a_start_exp + 1,
        a_end_exp - genomic,
    )
    valid = (pos >= 1) & (pos <= n_pos)
    if not valid.any():
        return matrix

    cols = pos[valid].astype(np.int64, copy=False) - 1
    rows_exp = rows_exp[valid]
    if use_scores:
        scores_exp = np.repeat(scores, lengths)[valid]
        np.add.at(matrix, (rows_exp, cols), scores_exp)
    else:
        np.add.at(matrix, (rows_exp, cols), 1)
    return matrix


def _bed_has_content(path: str) -> bool:
    """Cheap non-empty check that avoids loading the file."""
    import os
    try:
        return os.path.getsize(path) > 0
    except OSError:
        return False


def detect_nontrivial_bed_scores(
    xl_bed: str, n_sample: int = 5000
) -> bool:
    """Peek at the first ``n_sample`` rows of ``xl_bed`` and decide
    whether BED column 5 carries meaningful scores.

    Returns ``True`` if the score column contains more than one
    distinct value among the sampled rows AND at least one of those
    values is not 0 or 1 (i.e. it is not just a "presence" indicator).
    Used by the pipeline to nudge users toward ``--xl_score`` when
    they pass a track with real scores but leave the flag at its
    backward-compatible default of ``ignore``.

    A non-fatal sniffer: bad columns / weird formats just return
    ``False`` and the caller falls back to the legacy count-based
    behaviour.
    """
    try:
        df = pd.read_csv(
            xl_bed, sep='\t', header=None, comment='#',
            usecols=[4], names=['score'],
            dtype={'score': str}, nrows=n_sample,
        )
    except (ValueError, pd.errors.ParserError, OSError):
        return False
    if df.empty:
        return False
    # Try to coerce to numeric; non-numeric scores are ignored.
    nums = pd.to_numeric(df['score'], errors='coerce').dropna().to_numpy()
    if nums.size == 0:
        return False
    distinct = np.unique(nums)
    if distinct.size <= 1:
        return False
    non_indicator = ~np.isin(distinct, [0.0, 1.0])
    return bool(non_indicator.any())


def aggregate_legacy_columns(
    region,
    exon_categories: pd.Series,
    control_label: str = "control",
) -> pd.DataFrame:
    """Produce the legacy per-(category, position) aggregate frame.

    Used by ``bootstrap_contrast``, ``permutation_z`` and ``cluster_perm``
    to attach ``coverage`` / ``norm_coverage`` / ``fold_change`` columns
    to their plot frames so downstream TSV / plot helpers keep working.

    Accepts either a :class:`RegionCoverage` (fast, matrix-based) or the
    legacy long-form per-exon DataFrame (kept for backward compat with
    existing tests).
    """
    if isinstance(region, pd.DataFrame):
        agg = region.groupby(['name', 'position'], as_index=False).agg(
            coverage=('coverage', 'sum')
        )
        counts = pd.DataFrame({
            'name': exon_categories.index,
            'number_exons': exon_categories.values,
        })
        agg = agg.merge(counts, on='name', how='left')
    else:
        sums = region.category_sums()
        n_per_cat = {
            cat: int((region.exon_names == cat).sum()) for cat in sums
        }
        rows = []
        for cat, vec in sums.items():
            rows.append(pd.DataFrame({
                'name': cat,
                'position': region.positions,
                'coverage': vec.astype(np.int64),
                'number_exons': n_per_cat[cat],
            }))
        agg = (
            pd.concat(rows, ignore_index=True) if rows
            else pd.DataFrame(columns=['name', 'position', 'coverage',
                                       'number_exons'])
        )

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


def _fisher_linegraph_from_matrix(
    region: RegionCoverage,
    smoothing: int,
    control_label: str = "control",
) -> pd.DataFrame:
    """Reproduce the legacy per-(category, position) Fisher linegraph.

    Mirrors the schema the old long-form path emitted so the legacy
    ``--enrichment fisher`` plot, the TSV outputs, and the
    ``_aggregate_legacy_columns`` helpers used by other enrichment
    methods all keep working with no schema change.
    """
    sums = region.category_sums()
    counts = {
        cat: int((region.exon_names == cat).sum())
        for cat in sums
    }

    ctrl_sum = sums.get(control_label,
                        np.zeros(region.n_positions, dtype=np.int64))
    ctrl_n = counts.get(control_label, 0)
    ctrl_norm = np.zeros(region.n_positions, dtype=np.float64)
    if ctrl_n > 0:
        ctrl_norm = ctrl_sum.astype(np.float64) / ctrl_n
    ctrl_norm_safe = np.where(ctrl_norm == 0, 1e-6, ctrl_norm)

    rows = []
    for cat, cov in sums.items():
        n_exons = counts[cat]
        norm = (cov.astype(np.float64) / n_exons) if n_exons else np.zeros_like(
            cov, dtype=np.float64
        )
        fold = norm / ctrl_norm_safe

        # Per-position 2x2 Fisher: (cat hit, cat miss, ctrl hit, ctrl miss).
        cat_miss = np.maximum(n_exons - cov, 0)
        ctrl_miss = np.maximum(ctrl_n - ctrl_sum, 0)
        pvals = np.empty(region.n_positions, dtype=np.float64)
        # Tight Python loop over positions only (601 calls, not 49M rows).
        for i in range(region.n_positions):
            table = np.array([
                [cov[i], cat_miss[i]],
                [ctrl_sum[i], ctrl_miss[i]],
            ])
            pvals[i] = stats.fisher_exact(table)[1]
        # Signed -log10(p) matches the prior behaviour exactly.
        log10p = np.log10(1.0 / np.maximum(pvals, np.finfo(float).tiny))
        log10p = np.where(fold < 1, -log10p, log10p)
        smoothed = pd.Series(log10p).rolling(
            smoothing, center=True, win_type='gaussian'
        ).mean(std=2).to_numpy() if smoothing and smoothing > 1 else log10p.copy()

        rows.append(pd.DataFrame({
            'name': cat,
            'position': region.positions,
            'coverage': cov.astype(np.int64),
            'number_exons': n_exons,
            'norm_coverage': norm,
            'control_coverage': ctrl_sum.astype(np.int64),
            'control_number_exons': ctrl_n,
            'control_norm_coverage': ctrl_norm_safe,
            'fold_change': fold,
            'pvalue': pvals,
            '-log10pvalue': log10p,
            'label': region.label,
            '-log10pvalue_smoothed': smoothed,
        }))
    if not rows:
        return pd.DataFrame(columns=[
            'name', 'position', 'coverage', 'number_exons', 'norm_coverage',
            'control_coverage', 'control_number_exons',
            'control_norm_coverage', 'fold_change',
            'pvalue', '-log10pvalue', 'label', '-log10pvalue_smoothed',
        ])
    return pd.concat(rows, ignore_index=True)


def get_coverage_plot(
    xl_bed: str,
    df: pd.DataFrame,
    fai: str,
    window: int,
    exon_categories: pd.Series,
    label: str,
    smoothing: int = 15,
    score_mode: str = "ignore",
):
    """Compute per-exon CLIP coverage around one splice-site region.

    Parameters
    ----------
    xl_bed : str
        Path to (sorted) crosslink BED file.
    df : pd.DataFrame
        Splice-site BED frame produced by ``get_ss_bed`` (one row per
        exon-strand). ``name`` is ``category_chr:start-end;strand``.
    fai : str
        Path to the genome ``.fai`` (for ``slop -g``).
    window : int
        Half-window length in bases. Each exon gets a window of
        ``2 * window + 1`` positions.
    exon_categories : pd.Series
        Counts per category (for backward-compat schema only).
    label : str
        Splice-site region label.
    smoothing : int
        Gaussian rolling-mean window applied to the legacy Fisher
        ``-log10pvalue`` linegraph. Per-exon outputs are unsmoothed.
    score_mode : str
        How to use BED column 5 of ``xl_bed``. ``"ignore"`` (default)
        counts each overlap as ``+1``, matching the legacy behaviour.
        ``"raw"`` uses BED column 5 as-is; the per-transcript modes
        renormalise scores within each B.name (transcript) before
        accumulation. See :func:`_overlap_pairs_to_matrix` for details.

    Returns
    -------
    df_plot : pd.DataFrame
        Aggregated (category × position) frame with the legacy Fisher
        columns (``coverage``, ``norm_coverage``, ``fold_change``,
        ``pvalue``, ``-log10pvalue``, ``-log10pvalue_smoothed``, ...).
    region : RegionCoverage
        Compact per-exon coverage matrix and metadata for the
        per-exon enrichment methods (bootstrap_contrast,
        permutation_z, cluster_perm) and the heatmap.
    """
    ss_df, ss_pbt = _ss_df_to_pbt(df, window, fai)

    exon_names_full = ss_df['name'].to_numpy(dtype=object)
    name_to_row = {n: i for i, n in enumerate(exon_names_full)}
    n_exons = exon_names_full.size
    n_pos = 2 * window + 1

    # Split category vs exon-id once.
    split = np.array([n.split('_', 1) for n in exon_names_full], dtype=object)
    categories = split[:, 0].astype(object) if split.size else np.array(
        [], dtype=object
    )
    exon_ids = split[:, 1].astype(object) if split.size else np.array(
        [], dtype=object
    )

    xl_pbt = pbt.BedTool(xl_bed).sort()

    # Intersect: one row per (exon-window, crosslink) overlap pair. For
    # 380K crosslinks and ~80K exons the output is on the order of a
    # few million lines -- 50-100x smaller than ``coverage -d`` would
    # emit (which is one row per (exon, base) regardless of coverage).
    overlap = ss_pbt.intersect(
        xl_pbt, wa=True, wb=True, s=True, sorted=True, nonamecheck=True,
    )
    matrix = _overlap_pairs_to_matrix(
        overlap, name_to_row, n_exons, n_pos, score_mode=score_mode,
    )

    # Sort rows by (category, exon_id) so bootstrap / permutation
    # samples consume rows in a deterministic, RNG-stable order that
    # matches the legacy ``pivot_table(...).sort_index()`` behaviour
    # and produces reproducible CIs across runs.
    if n_exons:
        order = np.lexsort((exon_ids, categories))
        matrix = matrix[order]
        exon_ids = exon_ids[order]
        categories = categories[order]

    region = RegionCoverage(
        label=label,
        matrix=matrix,
        exon_ids=exon_ids,
        exon_names=categories,
        positions=np.arange(1, n_pos + 1, dtype=np.int32),
    )

    df_plot = _fisher_linegraph_from_matrix(region, smoothing)
    return df_plot, region

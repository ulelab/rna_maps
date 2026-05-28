"""Tests for BED-score-aware matrix building in rnamaps.coverage.

Exercises ``_overlap_pairs_to_matrix`` and ``detect_nontrivial_bed_scores``
without spinning up bedtools; the bedtools path is the same shape and is
covered indirectly by the end-to-end smoke test.
"""

from __future__ import annotations

import io
import os
import tempfile
from types import SimpleNamespace

import numpy as np

from rnamaps.coverage import (
    _overlap_pairs_to_matrix,
    detect_nontrivial_bed_scores,
)


def _write_overlap_lines(lines):
    """Persist a fake bedtools-intersect-wa-wb output to a temp file.

    Returns a stand-in object exposing the same ``.fn`` attribute that
    pybedtools' BedTool does.
    """
    tmp = tempfile.NamedTemporaryFile(
        mode='w', suffix='.bed', delete=False
    )
    for line in lines:
        tmp.write('\t'.join(str(c) for c in line) + '\n')
    tmp.close()
    return SimpleNamespace(fn=tmp.name)


def test_ignore_mode_counts_each_overlap_as_one():
    """Legacy behaviour: matrix is an int count of overlaps."""
    name_to_row = {"enhanced_e0": 0, "enhanced_e1": 1}
    n_exons = 2
    n_pos = 5

    # Two overlaps on row 0 at the same position; one on row 1.
    # A: chr1 start end name score strand  (col 0..5)
    # B: chr1 start end name score strand  (col 6..11)
    lines = [
        ['chr1', 100, 105, 'enhanced_e0', 0, '+',
         'chr1', 101, 102, 'b0', 0.7, '+'],
        ['chr1', 100, 105, 'enhanced_e0', 0, '+',
         'chr1', 101, 102, 'b1', 0.3, '+'],
        ['chr1', 200, 205, 'enhanced_e1', 0, '+',
         'chr1', 203, 204, 'b2', 5.0, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, n_exons, n_pos, score_mode='ignore'
        )
        assert m.dtype == np.int32
        # position is 1-indexed (b_start - a_start + 1 = 2)
        assert m[0, 1] == 2
        # second exon at pos 4
        assert m[1, 3] == 1
        assert m.sum() == 3
    finally:
        os.unlink(overlap.fn)


def test_raw_mode_accumulates_bed_scores():
    """``raw`` should sum the BED score column at each overlap."""
    name_to_row = {"enhanced_e0": 0}
    lines = [
        ['chr1', 100, 105, 'enhanced_e0', 0, '+',
         'chr1', 101, 102, 'b0', 0.7, '+'],
        ['chr1', 100, 105, 'enhanced_e0', 0, '+',
         'chr1', 101, 102, 'b1', 0.3, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, 1, 5, score_mode='raw'
        )
        assert m.dtype == np.float64
        assert m[0, 1] == 1.0  # 0.7 + 0.3
    finally:
        os.unlink(overlap.fn)


def test_per_transcript_zscore_zero_mean_unit_var_per_tx():
    """A per-transcript z-score within a single transcript should have
    zero mean (with ddof=0) and well-defined values."""
    name_to_row = {"e0": 0, "e1": 1, "e2": 2}
    lines = [
        ['chr1', 100, 110, 'e0', 0, '+',
         'chr1', 100, 101, 'txA', 1.0, '+'],
        ['chr1', 100, 110, 'e1', 0, '+',
         'chr1', 100, 101, 'txA', 2.0, '+'],
        ['chr1', 100, 110, 'e2', 0, '+',
         'chr1', 100, 101, 'txA', 3.0, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, 3, 10,
            score_mode='per_transcript_zscore',
        )
        # Z-scores of [1, 2, 3] (ddof=0) are [-sqrt(1.5), 0, +sqrt(1.5)].
        expected = np.array([-1.0, 0.0, 1.0]) / np.std(
            [1.0, 2.0, 3.0], ddof=0
        )
        np.testing.assert_allclose(
            [m[0, 0], m[1, 0], m[2, 0]], expected, atol=1e-12
        )
    finally:
        os.unlink(overlap.fn)


def test_multibase_feature_spreads_across_window_plus_strand():
    """A multi-base B-feature should deposit a count at every base it
    covers in the A-window, not just at its leftmost base.

    Regression test for the bug where finemapped 75-nt windows were
    only depositing +1 at b_start.
    """
    name_to_row = {"enhanced_e0": 0}
    n_exons = 1
    n_pos = 20

    # A-window: chr1:100-120 (so positions 1..20 map to genomic 100..119)
    # B-feature: chr1:103-108 (5 bases, should hit transcript pos 4..8)
    lines = [
        ['chr1', 100, 120, 'enhanced_e0', 0, '+',
         'chr1', 103, 108, 'b0', 1, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, n_exons, n_pos, score_mode='ignore'
        )
        # 5 contiguous positions should each have count 1.
        expected = np.zeros(n_pos, dtype=np.int32)
        expected[3:8] = 1  # 1-indexed pos 4..8 → 0-indexed cols 3..7
        np.testing.assert_array_equal(m[0], expected)
        assert m.sum() == 5
    finally:
        os.unlink(overlap.fn)


def test_multibase_feature_spreads_across_window_minus_strand():
    """Same as the +strand test but on - strand: bases should map to
    the mirrored transcript positions so + and - strand exons share a
    coordinate system in the matrix.
    """
    name_to_row = {"silenced_e0": 0}
    n_exons = 1
    n_pos = 20

    # A-window: chr1:100-120, strand -. Transcript pos 1 = genomic 119.
    # B-feature: chr1:103-108 (genomic bases 103..107).
    # For - strand: pos = a_end - genomic, so genomic 107..103 → pos 13..17.
    lines = [
        ['chr1', 100, 120, 'silenced_e0', 0, '-',
         'chr1', 103, 108, 'b0', 1, '-'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, n_exons, n_pos, score_mode='ignore'
        )
        expected = np.zeros(n_pos, dtype=np.int32)
        expected[12:17] = 1  # 1-indexed pos 13..17 → cols 12..16
        np.testing.assert_array_equal(m[0], expected)
        assert m.sum() == 5
    finally:
        os.unlink(overlap.fn)


def test_multibase_feature_clipped_to_window():
    """A B-feature that hangs off the end of the A-window should have
    only its in-window bases counted. bedtools -wa -wb reports the
    original (unclipped) B coords, so the matrix builder has to clip.
    """
    name_to_row = {"enhanced_e0": 0}
    n_exons = 1
    n_pos = 10

    # A-window: chr1:100-110 (positions 1..10 = genomic 100..109)
    # B-feature: chr1:107-115 (extends past the window by 5 bases).
    # Only genomic 107..109 (3 bases) are in-window → transcript pos 8..10.
    lines = [
        ['chr1', 100, 110, 'enhanced_e0', 0, '+',
         'chr1', 107, 115, 'b0', 1, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, n_exons, n_pos, score_mode='ignore'
        )
        expected = np.zeros(n_pos, dtype=np.int32)
        expected[7:10] = 1
        np.testing.assert_array_equal(m[0], expected)
        assert m.sum() == 3
    finally:
        os.unlink(overlap.fn)


def test_multibase_feature_raw_score_replicated_across_bases():
    """In ``raw`` mode each base of the B-feature should receive the
    feature's BED score (replicated, not divided by length).
    """
    name_to_row = {"enhanced_e0": 0}
    n_exons = 1
    n_pos = 10

    # B-feature: 3 bases wide, score 2.5 → each base gets 2.5.
    lines = [
        ['chr1', 100, 110, 'enhanced_e0', 0, '+',
         'chr1', 102, 105, 'b0', 2.5, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, n_exons, n_pos, score_mode='raw'
        )
        expected = np.zeros(n_pos, dtype=np.float64)
        expected[2:5] = 2.5
        np.testing.assert_allclose(m[0], expected)
    finally:
        os.unlink(overlap.fn)


def test_detect_nontrivial_bed_scores_count_data_is_trivial():
    """All-ones BED score column should be flagged as trivial."""
    tmp = tempfile.NamedTemporaryFile(
        mode='w', suffix='.bed', delete=False
    )
    for i in range(50):
        tmp.write(f"chr1\t{1000 + i}\t{1001 + i}\t.\t1\t+\n")
    tmp.close()
    try:
        assert detect_nontrivial_bed_scores(tmp.name) is False
    finally:
        os.unlink(tmp.name)


def test_detect_nontrivial_bed_scores_continuous_is_flagged():
    """A continuous BED score column should be flagged as non-trivial."""
    tmp = tempfile.NamedTemporaryFile(
        mode='w', suffix='.bed', delete=False
    )
    for i in range(50):
        tmp.write(
            f"chr1\t{1000 + i}\t{1001 + i}\t.\t{0.01 * (i + 1):.4f}\t+\n"
        )
    tmp.close()
    try:
        assert detect_nontrivial_bed_scores(tmp.name) is True
    finally:
        os.unlink(tmp.name)

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


def test_per_transcript_sum1_renormalises_b_names():
    """Scores within a single B.name should renormalise to sum to 1
    across the overlap rows that share that name, before accumulating
    onto the exon matrix.
    """
    name_to_row = {"enhanced_e0": 0, "enhanced_e1": 1}
    n_exons = 2
    n_pos = 10

    # Two crosslinks share txA (b_name); they should each contribute
    # 5/(5+5) = 0.5 to their respective positions. Without renorm they
    # would each contribute 5.0.
    lines = [
        ['chr1', 100, 110, 'enhanced_e0', 0, '+',
         'chr1', 100, 101, 'txA', 5.0, '+'],
        ['chr1', 200, 210, 'enhanced_e1', 0, '+',
         'chr1', 200, 201, 'txA', 5.0, '+'],
    ]
    overlap = _write_overlap_lines(lines)
    try:
        m = _overlap_pairs_to_matrix(
            overlap, name_to_row, n_exons, n_pos,
            score_mode='per_transcript_sum1',
        )
        # pos = b_start - a_start + 1 = 1 for both rows
        assert m[0, 0] == 0.5
        assert m[1, 0] == 0.5
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

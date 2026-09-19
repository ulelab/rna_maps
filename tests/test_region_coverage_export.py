"""Regression tests for event-level region-matrix exports."""

from __future__ import annotations

import numpy as np
import pandas as pd

from rnamaps.coverage import RegionCoverage, _ss_df_to_pbt


def test_save_tsv_preserves_duplicate_exon_coordinate_events(tmp_path):
    """Duplicate coordinates must remain two rows, not form a merge product."""
    exon_id = "chr1:100-200;+"
    region = RegionCoverage(
        label="middle_3ss",
        matrix=np.array([[1, 0, 0], [0, 0, 2]], dtype=np.int32),
        exon_ids=np.array([exon_id, exon_id]),
        exon_names=np.array(["control", "control"]),
        positions=np.array([1, 2, 3]),
    )
    heatmap_order = pd.DataFrame({
        "category": ["control", "control"],
        "exon_id": [exon_id, exon_id],
        "heatmap_signal": [1.0, 2.0],
        "heatmap_row": pd.array([2, 1], dtype="Int64"),
    })

    path = tmp_path / "region.tsv"
    region.save_tsv(path, heatmap_order)
    exported = pd.read_csv(path, sep="\t")

    assert len(exported) == 2
    assert exported["event_occurrence"].tolist() == [2, 1]
    assert exported["heatmap_row"].tolist() == [1, 2]
    assert exported["heatmap_signal"].tolist() == [2.0, 1.0]
    assert exported["total_crosslinks"].tolist() == [2, 1]


def test_splice_site_bed_uses_unique_internal_names_for_duplicate_events(
    monkeypatch,
):
    """Coverage rows must not collapse when event coordinates are duplicated."""
    frame = pd.DataFrame({
        "chr": ["chr1", "chr1"],
        "start": [100, 100],
        "end": [101, 101],
        "name": ["control_chr1:100-200;+", "control_chr1:100-200;+"],
        "score": [0.001, 0.001],
        "strand": ["+", "+"],
    })
    captured = {}

    class FakeBedTool:
        def sort(self):
            return self

        def slop(self, **kwargs):
            captured["slop"] = kwargs
            return self

    def fake_from_dataframe(bed_frame):
        captured["frame"] = bed_frame.copy()
        return FakeBedTool()

    monkeypatch.setattr(
        "rnamaps.coverage.pbt.BedTool.from_dataframe", fake_from_dataframe
    )

    indexed, _ = _ss_df_to_pbt(frame, window=50, fai="genome.fai")

    assert indexed["name"].tolist() == frame["name"].tolist()
    assert indexed["_coverage_row_id"].tolist() == ["0", "1"]
    assert captured["frame"]["name"].tolist() == ["0", "1"]
    assert captured["slop"] == {
        "l": 50, "r": 50, "s": True, "g": "genome.fai"
    }

"""Input-contract tests for VastDB cassette-exon lists."""

from __future__ import annotations

import pandas as pd
import pytest

from rnamaps.io_vastdb import load_vastdb_data


def _write_ids(tmp_path, name, ids):
    path = tmp_path / name
    path.write_text("".join(f"{event}\n" for event in ids))
    return path


def _write_annotation(tmp_path, rows):
    path = tmp_path / "EVENT_INFO.tab"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _annotation_row(event, central, upstream, downstream):
    chrom, start, end, strand = central
    return {
        "EVENT": event,
        "GENE": "GENE1",
        "COORD_o": f"{chrom}:{start}-{end}",
        "REF_CO": f"{chrom}:{strand}",
        "CO_C1": f"{chrom}:{upstream[0]}-{upstream[1]}",
        "CO_C2": f"{chrom}:{downstream[0]}-{downstream[1]}",
    }


def test_rejects_non_exon_skipping_event_ids(tmp_path):
    enhanced = _write_ids(tmp_path, "enhanced.txt", ["HsaALTA0000001"])
    annotation = _write_annotation(tmp_path, [])

    with pytest.raises(ValueError, match="HsaEX events only"):
        load_vastdb_data(
            enhanced, None, None, None, annotation, chroms=[]
        )


def test_rejects_duplicate_event_id_across_categories(tmp_path):
    enhanced = _write_ids(tmp_path, "enhanced.txt", ["HsaEX0000001"])
    control = _write_ids(tmp_path, "control.txt", ["HsaEX0000001"])
    annotation = _write_annotation(tmp_path, [])

    with pytest.raises(ValueError, match="unique across all input lists"):
        load_vastdb_data(
            enhanced, None, control, None, annotation, chroms=[]
        )


def test_rejects_distinct_events_for_same_central_exon(tmp_path):
    enhanced = _write_ids(tmp_path, "enhanced.txt", ["HsaEX0000001"])
    control = _write_ids(tmp_path, "control.txt", ["HsaEX0000002"])
    annotation = _write_annotation(tmp_path, [
        _annotation_row(
            "HsaEX0000001",
            ("chr1", 101, 200, "+"),
            (1, 50),
            (250, 300),
        ),
        _annotation_row(
            "HsaEX0000002",
            ("chr1", 101, 200, "+"),
            (10, 60),
            (260, 310),
        ),
    ])

    with pytest.raises(ValueError, match="same central exon coordinate"):
        load_vastdb_data(
            enhanced, None, control, None, annotation, chroms=[]
        )


def test_accepts_unique_hsaex_events(tmp_path):
    enhanced = _write_ids(tmp_path, "enhanced.txt", ["HsaEX0000001"])
    control = _write_ids(tmp_path, "control.txt", ["HsaEX0000002"])
    annotation = _write_annotation(tmp_path, [
        _annotation_row(
            "HsaEX0000001",
            ("chr1", 101, 200, "+"),
            (1, 50),
            (250, 300),
        ),
        _annotation_row(
            "HsaEX0000002",
            ("chr1", 401, 500, "-"),
            (550, 600),
            (300, 350),
        ),
    ])

    loaded = load_vastdb_data(
        enhanced, None, control, None, annotation, chroms=[]
    )

    assert loaded["EVENT"].tolist() == ["HsaEX0000001", "HsaEX0000002"]
    assert loaded["exonStart_0base"].tolist() == [100, 400]
    assert loaded["category"].tolist() == ["enhanced", "control"]

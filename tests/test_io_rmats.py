"""Regression tests for rMATS loading, focused on per-exon deduplication."""

from __future__ import annotations

import pandas as pd

from rnamaps.io_rmats import load_rmats_data


# Thresholds matching the CLI defaults (-xc/-xm/-xi/-xf/-xe/-ms).
DEFAULTS = dict(
    min_ctrl=-0.05, max_ctrl=0.05, max_inclusion=0.9,
    max_fdr=0.1, max_enh=-0.05, min_sil=0.05,
    chroms=None, no_constitutive=False,
)


def _write_rmats(tmp_path, events):
    """Write a minimal rMATS SE table.

    ``events`` is a list of dicts with keys ``exon`` (an (start, end) pair),
    ``dPSI``, ``FDR`` and optionally ``flanks`` — separate flanking
    coordinates are what make two rows distinct events for the same exon.
    """
    rows = []
    for i, ev in enumerate(events):
        start, end = ev['exon']
        up, down = ev.get('flanks', (1000 + i * 10, 9000 + i * 10))
        rows.append({
            'ID': i,
            'GeneID': f'"ENSG{i:011d}"',
            'geneSymbol': f'"GENE{i}"',
            'chr': ev.get('chr', 'chr1'),
            'strand': ev.get('strand', '+'),
            'exonStart_0base': start,
            'exonEnd': end,
            'upstreamES': up,
            'upstreamEE': up + 100,
            'downstreamES': down,
            'downstreamEE': down + 100,
            'FDR': ev['FDR'],
            'IncLevel1': '0.5,0.5',
            'IncLevel2': '0.5,0.5',
            'IncLevelDifference': ev['dPSI'],
        })
    path = tmp_path / 'rmats.txt'
    pd.DataFrame(rows).to_csv(path, sep='\t', index=False)
    return str(path)


def test_dedup_keeps_single_row_for_two_way_tie(tmp_path):
    """A 2-way tie at max |dPSI| must not leak two rows for one exon.

    The averaged rank of a 2-way tie is 1.5, so the old ``rank < 2`` filter
    kept both rows.
    """
    path = _write_rmats(tmp_path, [
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01},
        {'exon': (100, 200), 'dPSI': -0.2, 'FDR': 0.01},
        {'exon': (100, 200), 'dPSI': 0.01, 'FDR': 1.0},
    ])

    df = load_rmats_data(path, **DEFAULTS)

    assert len(df) == 1
    # Ties break on original row order, so the first of the two wins.
    assert df.iloc[0]['dPSI'] == 0.2
    assert df.iloc[0]['category'] == 'silenced'


def test_dedup_keeps_exon_with_three_way_tie(tmp_path):
    """A >=3-way tie ranked 2.0 and dropped the exon entirely."""
    path = _write_rmats(tmp_path, [
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01},
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01},
        {'exon': (100, 200), 'dPSI': -0.2, 'FDR': 0.01},
    ])

    df = load_rmats_data(path, **DEFAULTS)

    assert len(df) == 1
    assert df.iloc[0]['exonStart_0base'] == 100


def test_dedup_keeps_most_extreme_dpsi(tmp_path):
    """The retained event is still the most extreme one, sign-agnostic."""
    path = _write_rmats(tmp_path, [
        {'exon': (100, 200), 'dPSI': 0.02, 'FDR': 1.0},
        {'exon': (100, 200), 'dPSI': -0.30, 'FDR': 0.01},
        {'exon': (100, 200), 'dPSI': 0.10, 'FDR': 0.01},
    ])

    df = load_rmats_data(path, **DEFAULTS)

    assert len(df) == 1
    assert df.iloc[0]['dPSI'] == -0.30
    assert df.iloc[0]['category'] == 'enhanced'


def test_distinct_exons_are_all_retained_in_input_order(tmp_path):
    """Deduplication must not reorder or drop genuinely distinct exons."""
    path = _write_rmats(tmp_path, [
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01},
        {'exon': (300, 400), 'dPSI': 0.0, 'FDR': 1.0},
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01},
        {'exon': (500, 600), 'dPSI': -0.2, 'FDR': 0.01},
    ])

    df = load_rmats_data(path, **DEFAULTS)

    assert df['exonStart_0base'].tolist() == [100, 300, 500]
    assert df['category'].tolist() == ['silenced', 'control', 'enhanced']


def test_same_coordinates_on_opposite_strands_are_separate_exons(tmp_path):
    """The exon key includes strand, so both strands survive."""
    path = _write_rmats(tmp_path, [
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01, 'strand': '+'},
        {'exon': (100, 200), 'dPSI': 0.2, 'FDR': 0.01, 'strand': '-'},
    ])

    df = load_rmats_data(path, **DEFAULTS)

    assert sorted(df['strand'].tolist()) == ['+', '-']

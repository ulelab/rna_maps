"""Track 1 input loader: rMATS differential splicing output."""

import logging

import numpy as np
import pandas as pd


def load_rmats_data(de_file, min_ctrl, max_ctrl, max_inclusion,
                    max_fdr, max_enh, min_sil, chroms, no_constitutive):
    """
    Load and categorise exons from rMATS output.

    Reads rMATS SE file, computes maxPSI, deduplicates, assigns categories
    from dPSI/FDR thresholds, and corrects upstream/downstream labels for
    minus-strand genes (rMATS labels by genomic position, not transcript order).

    Returns DataFrame with canonical columns for the shared pipeline.
    """
    logging.info("=" * 60)
    logging.info("INPUT MODE: rMATS")
    logging.info("=" * 60)

    rmats = pd.read_csv(de_file, sep='\t')

    if 'exonStart_0base' not in rmats.columns:
        raise ValueError(
            "Input file does not appear to be rMATS format "
            "(missing 'exonStart_0base' column)"
        )

    rmats = rmats[rmats['chr'].isin(chroms)]

    # Compute max PSI across all samples
    rmats['inclusion'] = (
        rmats.IncLevel1.str.split(',') + rmats.IncLevel2.str.split(',')
    )
    rmats['inclusion'] = rmats['inclusion'].apply(
        lambda x: max([float(y) for y in x if y != 'NA'])
    )

    df_rmats = rmats.loc[:, [
        'chr', 'exonStart_0base', 'exonEnd', 'FDR', 'IncLevelDifference',
        'strand', 'inclusion',
        'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'
    ]].rename(columns={
        'IncLevelDifference': 'dPSI',
        'inclusion': 'maxPSI'
    }).reset_index()

    # Deduplicate: keep the most extreme dPSI per exon
    mask = df_rmats.groupby(
        ['chr', 'exonStart_0base', 'exonEnd', 'strand']
    )['dPSI'].transform(lambda x: abs(x).rank(ascending=False)) < 2
    df_rmats = df_rmats[mask]

    # Assign categories from thresholds
    conditions = [
        (df_rmats["dPSI"].gt(min_sil) & df_rmats["FDR"].lt(max_fdr)),      # silenced
        (df_rmats["dPSI"].lt(max_enh) & df_rmats["FDR"].lt(max_fdr)),      # enhanced
        (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl)
         & df_rmats["maxPSI"].gt(max_inclusion)),                            # constitutive
        (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl)),   # control
    ]
    choices = ["silenced", "enhanced", "constituitive", "control"]
    df_rmats["category"] = np.select(conditions, choices, default=None)

    # Filter out constitutive if requested
    if no_constitutive:
        df_rmats = df_rmats[df_rmats['category'] != 'constituitive']

    # ---------------------------------------------------------------
    # FIX: Swap upstream/downstream for minus-strand genes.
    #
    # rMATS labels "upstream" and "downstream" by genomic coordinate
    # (lower = upstream, higher = downstream). For minus-strand genes
    # this is inverted relative to transcript order. We swap here so
    # that after loading, upstream/downstream always mean transcript-
    # relative direction, matching VastDB's CO_C1/CO_C2 convention.
    # ---------------------------------------------------------------
    minus_mask = df_rmats['strand'] == '-'
    cols_to_swap = [('upstreamES', 'downstreamES'), ('upstreamEE', 'downstreamEE')]
    for col_a, col_b in cols_to_swap:
        df_rmats.loc[minus_mask, [col_a, col_b]] = \
            df_rmats.loc[minus_mask, [col_b, col_a]].values

    logging.info(f"Loaded {len(df_rmats)} categorised exons from rMATS")
    logging.info(f"Category distribution:\n{df_rmats['category'].value_counts()}")

    return df_rmats

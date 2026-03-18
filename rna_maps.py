#!/usr/bin/env python3
"""
RNA Maps - Unified Script

Generates RNA maps showing RBP binding patterns around regulated exons.
Supports two input modes:
  1. rMATS mode:  Takes rMATS differential splicing output, auto-categorises exons
  2. VastDB mode: Takes pre-curated VastDB EVENT ID lists + EVENT_INFO annotation

Key features:
- 6 splice site regions (upstream_3ss, upstream_5ss, middle_3ss, middle_5ss,
  downstream_3ss, downstream_5ss)
- Coverage plots with Fisher's exact test enrichment
- Per-exon heatmaps
- Exon length distributions
- Optional multivalency analysis (rMATS mode, requires germs.R)

Usage examples:

  # rMATS mode
  python rna_maps.py -i rMATS_SE.txt -x CLIP.bed -f hg19.fa -fi hg19.fa.fai -o output -p PTBP1

  # VastDB mode
  python rna_maps.py --vastdb_mode \\
    --vastdb_enhanced enhanced_ids.txt \\
    --vastdb_silenced silenced_ids.txt \\
    --vastdb_control control_ids.txt \\
    --vastdb_constitutive constitutive_ids.txt \\
    --vastdb_annotation EVENT_INFO-hg38.tab \\
    -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p AQR_K562
"""

import matplotlib
import matplotlib.ticker as mticker
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from scipy.ndimage import gaussian_filter1d
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import pybedtools as pbt
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import os
import argparse
import random
import string
import logging
import datetime
import time
import re


# ===============================================================================
# CONFIGURATION
# ===============================================================================

sns.set_style("whitegrid", {'legend.frameon': True})

colors_dict = {
    'all': '#D3D3D3',
    'ctrl': '#408F76',
    'enh': '#F30C08',
    'sil': '#005299',
    'enhrest': '#FFB122',
    'silrest': '#6DC2F5',
    'const': '#666666'
}
linewidth = 3
dashes = False


# ===============================================================================
# UTILITY FUNCTIONS
# ===============================================================================

def setup_logging(output_path):
    """Sets up logging to file and console."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"execution_{timestamp}.log"
    log_path = os.path.join(output_path, log_filename)

    logger = logging.getLogger()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    file_handler.flush = file_handler.stream.flush

    start_time = time.time()

    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
    logging.getLogger('matplotlib.backends').setLevel(logging.ERROR)
    logging.getLogger('fontTools').setLevel(logging.ERROR)
    logging.getLogger('fontTools.subset').setLevel(logging.ERROR)

    logger.info(f"Starting script execution at {timestamp}")

    return log_filename, start_time, logger


def log_runtime(start_time, logger):
    """Logs the script execution time."""
    end_time = time.time()
    total_seconds = int(end_time - start_time)
    minutes, seconds = divmod(total_seconds, 60)
    runtime_str = f"{minutes}m {seconds}s" if minutes > 0 else f"{seconds}s"
    logger.info(f"Script completed in {runtime_str}.")
    for handler in logger.handlers:
        handler.flush()


def smooth_coverage(df, window_size=10, std=2):
    """Smooth coverage data using rolling Gaussian window."""
    result = df.copy()
    groups = []

    for (exon_id, label), group_df in result.groupby(['exon_id', 'label']):
        if len(group_df) < 3:
            groups.append(group_df)
            continue

        group_sorted = group_df.sort_values('position')

        if len(group_sorted) >= window_size:
            values = group_sorted['coverage'].values
            s = pd.Series(values)

            smoothed = s.rolling(
                window=window_size,
                center=True,
                win_type='gaussian'
            ).mean(std=std)

            smoothed = smoothed.fillna(s)
            group_sorted['coverage'] = smoothed.values

        groups.append(group_sorted)

    return pd.concat(groups)


# ===============================================================================
# TRACK 1: rMATS DATA LOADING
# ===============================================================================

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


# ===============================================================================
# TRACK 2: VASTDB DATA LOADING
# ===============================================================================

def parse_event_info_coordinates(event_info_df):
    """
    Parse genomic coordinates from EVENT_INFO file.
    Returns DataFrame with exon and flanking exon coordinates.
    """
    coords_list = []

    for idx, row in event_info_df.iterrows():
        coord_str = str(row['COORD_o'])
        match = re.match(r'(chr[\w]+):(\d+)-(\d+)', coord_str)
        if not match:
            continue

        ref_co = str(row.get('REF_CO', ''))
        strand_match = re.search(r':([+-])$', ref_co)
        strand = strand_match.group(1) if strand_match else '+'

        upstream_start, upstream_end = np.nan, np.nan
        downstream_start, downstream_end = np.nan, np.nan

        if pd.notna(row.get('CO_C1', np.nan)):
            co_c1_match = re.match(r'chr[\w]+:(\d+)-(\d+)', str(row['CO_C1']))
            if co_c1_match:
                upstream_start = int(co_c1_match.group(1))
                upstream_end = int(co_c1_match.group(2))

        if pd.notna(row.get('CO_C2', np.nan)):
            co_c2_match = re.match(r'chr[\w]+:(\d+)-(\d+)', str(row['CO_C2']))
            if co_c2_match:
                downstream_start = int(co_c2_match.group(1))
                downstream_end = int(co_c2_match.group(2))

        coords_list.append({
            'EVENT': row['EVENT'],
            'GENE': row['GENE'],
            'chr': match.group(1),
            'exonStart_1based': int(match.group(2)),
            'exonEnd': int(match.group(3)),
            'strand': strand,
            'upstreamES': upstream_start,
            'upstreamEE': upstream_end,
            'downstreamES': downstream_start,
            'downstreamEE': downstream_end,
        })

    coords_df = pd.DataFrame(coords_list)

    if len(coords_df) > 0:
        logging.info(f"Parsed coordinates for {len(coords_df)} events")

    return coords_df


def load_vastdb_data(enhanced_file, silenced_file, control_file,
                     constitutive_file, event_info_file, chroms):
    """
    Load VastDB ID lists and look up coordinates from EVENT_INFO.

    Categories are assigned by list membership. Coordinates are converted
    from 1-based (VastDB) to 0-based (BED). CO_C1/CO_C2 are already in
    transcript order, so no strand correction is needed.

    Returns DataFrame with canonical columns for the shared pipeline.
    """
    logging.info("=" * 60)
    logging.info("INPUT MODE: VastDB ID Lists")
    logging.info("=" * 60)

    # Read ID lists
    def read_id_list(filepath, category):
        if filepath is None:
            return []
        with open(filepath) as f:
            ids = [line.strip() for line in f
                   if line.strip() and not line.startswith('#')]
        logging.info(f"Loaded {len(ids)} {category} IDs")
        return ids

    enhanced_ids = read_id_list(enhanced_file, 'enhanced')
    silenced_ids = read_id_list(silenced_file, 'silenced')
    control_ids = read_id_list(control_file, 'control')
    constitutive_ids = read_id_list(constitutive_file, 'constitutive')

    all_ids = enhanced_ids + silenced_ids + control_ids + constitutive_ids

    if len(all_ids) == 0:
        raise ValueError("No EVENT IDs provided in any file!")

    logging.info(f"Total EVENT IDs: {len(all_ids)}")

    # Create ID to category mapping
    id_to_category = {}
    for eid in enhanced_ids:
        id_to_category[eid] = 'enhanced'
    for eid in silenced_ids:
        id_to_category[eid] = 'silenced'
    for eid in control_ids:
        id_to_category[eid] = 'control'
    for eid in constitutive_ids:
        id_to_category[eid] = 'constituitive'  # Match original spelling

    # Load EVENT_INFO
    logging.info(f"Loading EVENT_INFO: {event_info_file}")
    event_info = pd.read_csv(event_info_file, sep='\t')
    logging.info(f"EVENT_INFO contains {len(event_info)} events")

    # Parse coordinates
    event_coords = parse_event_info_coordinates(event_info)

    # Filter to our IDs
    event_coords_subset = event_coords[event_coords['EVENT'].isin(all_ids)].copy()

    n_matched = len(event_coords_subset)
    n_total = len(all_ids)
    logging.info(f"Matched {n_matched}/{n_total} IDs to coordinates "
                 f"({100 * n_matched / n_total:.1f}%)")

    if n_matched == 0:
        raise ValueError("No EVENT IDs matched to coordinates!")

    # Convert 1-based to 0-based
    event_coords_subset['exonStart_0base'] = event_coords_subset['exonStart_1based'] - 1
    event_coords_subset['upstreamES'] = event_coords_subset['upstreamES'] - 1
    event_coords_subset['downstreamES'] = event_coords_subset['downstreamES'] - 1

    # Assign categories
    event_coords_subset['category'] = event_coords_subset['EVENT'].map(id_to_category)

    # Add FDR placeholder (not used but needed for column contract)
    event_coords_subset['FDR'] = 0.001

    # Filter to valid chromosomes
    event_coords_subset = event_coords_subset[event_coords_subset['chr'].isin(chroms)]

    logging.info(f"Category distribution:\n{event_coords_subset['category'].value_counts()}")

    return event_coords_subset


# ===============================================================================
# SHARED PIPELINE FUNCTIONS
# ===============================================================================

def apply_subsetting(df_rmats, no_constitutive):
    """
    Subset control and constitutive exons to match the largest regulated
    category count. Returns (subsetted df, original_counts dict).
    """
    category_counts = df_rmats['category'].value_counts()
    original_counts = {cat: count for cat, count in category_counts.items()}

    target_count = 0
    if 'enhanced' in category_counts and 'silenced' in category_counts:
        target_count = max(category_counts['enhanced'], category_counts['silenced'])
    elif 'enhanced' in category_counts:
        target_count = category_counts['enhanced']
    elif 'silenced' in category_counts:
        target_count = category_counts['silenced']

    # Subset control
    if 'control' in category_counts and category_counts['control'] > target_count > 0:
        control_indices = df_rmats[df_rmats['category'] == 'control'].index
        control_keep = np.random.choice(control_indices, target_count, replace=False)
        drop_mask = df_rmats.index.isin(control_indices) & ~df_rmats.index.isin(control_keep)
        df_rmats = df_rmats[~drop_mask]
        logging.info(f"Subsetted control exons from {category_counts['control']} to {target_count}")

    # Subset constitutive
    if (not no_constitutive
            and 'constituitive' in category_counts
            and category_counts['constituitive'] > target_count > 0):
        const_indices = df_rmats[df_rmats['category'] == 'constituitive'].index
        const_keep = np.random.choice(const_indices, target_count, replace=False)
        drop_mask = df_rmats.index.isin(const_indices) & ~df_rmats.index.isin(const_keep)
        df_rmats = df_rmats[~drop_mask]
        logging.info(f"Subsetted constitutive exons from "
                     f"{category_counts['constituitive']} to {target_count}")

    return df_rmats, original_counts


def get_ss_bed(df, pos_col, neg_col):
    """
    Create BED file for splice sites (handles strand orientation).

    pos_col: column to use for + strand (one end of the exon)
    neg_col: column to use for - strand (other end of same exon)

    Both columns should reference the SAME exon. The upstream/downstream
    strand correction happens at load time, not here.
    """
    df = df.copy()

    df['exon_id'] = (
        df['category'] + "_" +
        df['chr'].astype(str) + ":" +
        df['exonStart_0base'].astype(str) + "-" +
        df['exonEnd'].astype(str) + ";" +
        df['strand'].astype(str)
    )

    ss_pos = df.loc[df['strand'] == "+",
                    ['chr', pos_col, pos_col, 'exon_id', 'FDR', 'strand']]
    ss_pos.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_pos.start = ss_pos.start.transform(lambda x: x - 1)

    ss_n = df.loc[df['strand'] == "-",
                  ['chr', neg_col, neg_col, 'exon_id', 'FDR', 'strand']]
    ss_n.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_n.end = ss_n.end.transform(lambda x: x + 1)

    ss = pd.concat([ss_pos, ss_n])

    return ss


def get_coverage_plot(xl_bed, df, fai, window, exon_categories, label,
                      smoothing=15):
    """
    Calculate coverage of CLIP crosslinks around splice sites with
    Fisher's exact test enrichment against control exons.
    """
    df = df.loc[df.name != "."]
    xl_bed = pbt.BedTool(xl_bed).sort()
    pbt_df = pbt.BedTool.from_dataframe(
        df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    ).sort().slop(l=window, r=window, s=True, g=fai)

    df_coverage = pbt_df.coverage(
        b=xl_bed, **{'sorted': True, 's': True, 'd': True, 'nonamecheck': True}
    ).to_dataframe()[['thickStart', 'thickEnd', 'strand', 'name']]

    df_coverage.rename(
        columns={'thickStart': 'position', 'thickEnd': 'coverage'}, inplace=True
    )

    df_plot = df_coverage

    # Adjust positions based on strand
    df_plot.loc[df_plot.strand == '+', 'position'] = \
        df_plot['position'].astype('int32')
    df_plot.loc[df_plot.strand == '-', 'position'] = \
        abs(2 * window + 2 - df_plot['position'])

    df_plot = df_plot.loc[df_plot.name != "."]

    df_plot = df_plot.join(
        df_plot.pop('name').str.split('_', n=1, expand=True).rename(
            columns={0: 'name', 1: 'exon_id'}
        )
    )

    heatmap_plot = df_plot.copy()
    heatmap_plot['label'] = label

    # Aggregate coverage
    df_plot = df_plot.groupby(
        ['name', 'position'], as_index=False
    ).agg({'coverage': 'sum'})

    exon_cat = pd.DataFrame({
        'name': exon_categories.index,
        'number_exons': exon_categories.values
    })
    df_plot = df_plot.merge(exon_cat, how="left")

    df_plot['norm_coverage'] = np.where(
        df_plot['coverage'] == 0,
        df_plot['coverage'],
        df_plot['coverage'] / df_plot['number_exons']
    )

    # Fisher test against control
    df_plot_ctrl = df_plot.loc[
        df_plot.name == "control"
    ][["position", "coverage", "number_exons"]]
    df_plot_ctrl.columns = ["position", "control_coverage", "control_number_exons"]
    df_plot = df_plot.merge(df_plot_ctrl, how="left")
    df_plot['control_norm_coverage'] = (
        df_plot["control_coverage"] / df_plot["control_number_exons"]
    )

    # Pseudo-count for zero values
    df_plot.loc[
        df_plot['control_norm_coverage'] == 0, ['control_norm_coverage']
    ] = 0.000001
    df_plot['fold_change'] = (
        df_plot["norm_coverage"] / df_plot["control_norm_coverage"]
    )

    # Statistical test
    contingency_table = list(zip(
        df_plot['coverage'],
        df_plot['number_exons'] - df_plot['coverage'],
        df_plot['control_coverage'],
        df_plot['control_number_exons'] - df_plot['control_coverage']
    ))
    contingency_table = [np.array(table).reshape(2, 2)
                         for table in contingency_table]
    df_plot['pvalue'] = [stats.fisher_exact(table)[1]
                         for table in contingency_table]

    df_plot['-log10pvalue'] = np.log10(1 / df_plot['pvalue'])
    df_plot['label'] = label

    df_plot.loc[
        df_plot['fold_change'] < 1, ['-log10pvalue']
    ] = df_plot['-log10pvalue'] * -1
    df_plot['-log10pvalue_smoothed'] = df_plot['-log10pvalue'].rolling(
        smoothing, center=True, win_type="gaussian"
    ).mean(std=2)

    return df_plot, heatmap_plot


def set_legend_text(legend, exon_categories, original_counts=None):
    """Set legend text with optional subset information."""
    exon_cat = pd.DataFrame({
        'name': exon_categories.index,
        'number_exons': exon_categories.values
    })
    categories = exon_cat['name'].unique()
    legend.set_bbox_to_anchor([1.08, 0.75])
    legend.set_title("")

    legend_idx = 0
    for category in ['constituitive', 'control', 'enhanced', 'silenced']:
        if category in categories:
            count = exon_cat[exon_cat['name'] == category]['number_exons'].values[0]
            text = f"{category.capitalize()} ({count}"
            if original_counts and category in original_counts:
                orig_count = original_counts[category]
                if orig_count > count:
                    text += f", subset from {orig_count}"
            text += ")"
            legend.texts[legend_idx].set_text(text)
            legend_idx += 1


def add_enrichment_marker(fig, ax):
    """Add enrichment/depletion marker with arrows."""
    y_min, y_max = ax.get_ylim()
    marker_height = 0.4
    normalized_bottom = (0 - y_min) / (y_max - y_min) - marker_height / 2

    marker_ax = fig.add_axes(
        [0.94, normalized_bottom, 0.06, marker_height],
        label=f"enrichment_marker_{ax.get_title()}"
    )
    marker_ax.set_xticks([])
    marker_ax.set_yticks([])
    for spine in marker_ax.spines.values():
        spine.set_visible(False)

    marker_ax.text(
        0.5, 0.5, "↑\nenriched\n\n\n↓\ndepleted",
        ha='center', va='center', fontsize=11,
        transform=marker_ax.transAxes
    )
    return marker_ax


def get_multivalency_scores(df, fai, window, genome_fasta, output_dir,
                            name, type, germsdir):
    """Return multivalency scores around df features extended by windows."""
    df = df.loc[
        (df.name != ".") & (pd.notnull(df.name)) & (df.name != "None")
    ]
    df = df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    df.columns = ['chr', 'start', 'stop', 'name', 'score', 'strand']
    df['name'] = df['name'].apply(
        lambda x: (str(x) + "XX" +
                   ''.join(random.choice(string.ascii_lowercase) for _ in range(6)))
    )
    pbt_df = pbt.BedTool.from_dataframe(
        df[['chr', 'start', 'stop', 'name', 'score', 'strand']]
    ).sort().slop(l=2 * window, r=2 * window, s=True, g=fai)
    logging.info("Number of sites: " + str(len(pbt_df)))
    pbts = pbt.BedTool.filter(
        pbt_df, lambda x: len(x) == (4 * window) + 1
    ).saveas()
    logging.info(
        "Number of seqs after filtering off-chrom: " + str(len(pbts))
    )
    pbts.sequence(fi=genome_fasta, name=True).save_seqs(
        f'{output_dir}/{name}_{type}_temp.fa'
    )
    logging.info("Running germs to calculate multivalency scores...")

    os.system(
        "RScript --vanilla " + germsdir + "/germs.R -f "
        + f'{output_dir}/{name}_{type}_temp.fa' + " -w 100 -s 20"
    )
    os.system("gunzip -f *multivalency.tsv.gz")
    mdf = pd.read_csv(
        f'{output_dir}/{name}_{type}_temp_5_101_21.multivalency.tsv',
        sep='\t', header=0
    )
    os.system(f'rm {output_dir}/{name}_{type}_temp_5_101_21.multivalency.tsv')
    os.system(f'rm {output_dir}/{name}_{type}_temp.fa')
    mdf['position'] = np.tile(np.arange(0, 4 * window - 3), len(pbts))
    mdf[['exon_type', 'label', 'roname']] = mdf['sequence_name'].str.split(
        r'XX|_', expand=True
    )

    # Top contributing kmers
    filtered_mdf = mdf[mdf['exon_type'].isin(['enhanced', 'silenced'])]
    enhanced_df = filtered_mdf[filtered_mdf['exon_type'] == 'enhanced']
    silenced_df = filtered_mdf[filtered_mdf['exon_type'] == 'silenced']

    def calculate_top_kmers(df):
        kmer_scores = df.groupby('kmer')['smoothed_kmer_multivalency'].sum()
        top_kmers = kmer_scores.nlargest(5).index
        return df[df['kmer'].isin(top_kmers)]

    top_kmers_enhanced = calculate_top_kmers(enhanced_df)
    top_kmers_silenced = calculate_top_kmers(silenced_df)
    top_kmers_df = pd.concat([top_kmers_enhanced, top_kmers_silenced])

    mdf = mdf.groupby(
        ['exon_type', 'position'], as_index=False
    ).agg({'smoothed_kmer_multivalency': 'mean'}).reset_index()
    top_kmers_df = top_kmers_df.groupby(
        ['position', 'exon_type', 'kmer'], as_index=False
    ).agg({'smoothed_kmer_multivalency': 'mean'}).reset_index()

    mdf['type'] = type
    top_kmers_df['type'] = type

    return mdf, top_kmers_df


# ===============================================================================
# SHARED PLOTTING FUNCTIONS
# ===============================================================================

def plot_exon_lengths(df_rmats, output_dir, FILEname):
    """Generate exon length box plots."""
    df_rmats["regulated_exon_length"] = (
        df_rmats['exonEnd'] - df_rmats['exonStart_0base']
    )
    df_rmats["first_exon_length"] = (
        df_rmats['upstreamEE'] - df_rmats['upstreamES']
    )
    df_rmats["second_exon_length"] = (
        df_rmats['downstreamEE'] - df_rmats['downstreamES']
    )
    df_rmats.loc[df_rmats.strand == '+', 'upstream_exon_length'] = \
        df_rmats["first_exon_length"]
    df_rmats.loc[df_rmats.strand == '-', 'upstream_exon_length'] = \
        df_rmats["second_exon_length"]
    df_rmats.loc[df_rmats.strand == '+', 'downstream_exon_length'] = \
        df_rmats["second_exon_length"]
    df_rmats.loc[df_rmats.strand == '-', 'downstream_exon_length'] = \
        df_rmats["first_exon_length"]

    exon_length_df = df_rmats[[
        "regulated_exon_length", "upstream_exon_length",
        "downstream_exon_length", "category"
    ]]
    exon_length_df = exon_length_df.melt(
        id_vars=["category"], var_name="exon_type", value_name="exon_length"
    )

    palette_exon_len = [
        colors_dict['ctrl'], colors_dict['const'], colors_dict['enh'],
        colors_dict['enhrest'], colors_dict['sil'], colors_dict['silrest'],
        colors_dict['all']
    ]

    sns.set(rc={'figure.figsize': (15, 5)})
    sns.set_style("whitegrid")
    g = sns.catplot(
        data=exon_length_df, x='category', y='exon_length', col='exon_type',
        kind='box', col_wrap=3, showfliers=False,
        col_order=["upstream_exon_length", "regulated_exon_length",
                    "downstream_exon_length"],
        order=["control", "constituitive", "enhanced", "enhanced_rest",
               "silenced", "silenced_rest"],
        palette=palette_exon_len, hue='category', legend=False
    )
    titles = ["Upstream Exon", "Middle Exon", "Downstream exon"]
    for ax, title in zip(g.axes.flat, titles):
        ax.set_title(title)
    g.set_xticklabels(rotation=45)
    g.set(xlabel=None)
    g.axes[0].set_ylabel('Exon length (bp)')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{FILEname}_exon_length.pdf')
    pbt.helpers.cleanup()
    logging.info(f"Saved exon length plot to {output_dir}/{FILEname}_exon_length.pdf")


def plot_heatmap(heat_df, exon_categories, window, all_sites,
                 output_dir, FILEname):
    """Generate per-exon heatmap from binary coverage data."""
    # Total exons covered table
    grouped_heat_df = heat_df.groupby(
        ['exon_id', 'label', 'name']
    )['coverage'].sum().reset_index()
    filtered_heat_df = grouped_heat_df[grouped_heat_df['coverage'] > 0]
    count_heat_df = filtered_heat_df.groupby(
        ['label', 'name']
    )['exon_id'].nunique().reset_index()
    count_heat_df.rename(columns={'exon_id': 'exon_count'}, inplace=True)
    final_heat_df = count_heat_df.pivot(
        index='name', columns='label', values='exon_count'
    ).fillna(0).reset_index()
    exon_categories_df = exon_categories.reset_index()
    exon_categories_df.columns = ['name', 'total_exons_after_subsetting']
    final_heat_df = final_heat_df.merge(exon_categories_df, on='name', how='left')
    final_heat_df.to_csv(
        f'{output_dir}/{FILEname}_totalExonsCovered.tsv', sep="\t", index=False
    )

    # Binarise and smooth
    df = heat_df.copy()
    df['coverage'] = (df['coverage'] > 0).astype(int)
    df = smooth_coverage(df)

    if not all_sites:
        labels = ['upstream_5ss', 'middle_3ss', 'middle_5ss', 'downstream_3ss']
    else:
        labels = ['upstream_3ss', 'upstream_5ss', 'middle_3ss', 'middle_5ss',
                  'downstream_3ss', 'downstream_5ss']

    # Total signal per exon
    exon_totals = df.groupby('exon_id')['coverage'].sum().reset_index()
    exon_totals.columns = ['exon_id', 'total_signal']

    # Remove exons with no signal
    exons_with_signal = exon_totals[exon_totals['total_signal'] > 0]['exon_id']
    df = df[df['exon_id'].isin(exons_with_signal)]

    exon_names = df[['exon_id', 'name']].drop_duplicates(
        subset=['exon_id']
    ).set_index('exon_id')['name']

    label_data = {}
    for label in labels:
        label_df = df[df['label'] == label]
        if len(label_df) == 0:
            continue
        pivot = label_df.pivot_table(
            index='exon_id', columns='position',
            values='coverage', fill_value=0
        )
        label_data[label] = pivot

    # Common exons
    common_exons = set()
    first = True
    for label, pivot in label_data.items():
        if first:
            common_exons = set(pivot.index)
            first = False
        else:
            common_exons = common_exons.union(set(pivot.index))

    if len(common_exons) == 0:
        logging.info("No exons with signal for heatmap — skipping")
        return

    exon_info = exon_totals.set_index('exon_id').loc[list(common_exons)]
    exon_info['name'] = exon_names.loc[exon_info.index]
    exon_info = exon_info.sort_values(['name', 'total_signal'], ascending=[True, False])
    sorted_exon_ids = exon_info.index.tolist()

    # Set up figure
    width = max(15, len(labels) * 4)
    height = max(3, len(sorted_exon_ids) * 0.002)
    fig = plt.figure(figsize=(width, height))
    fig.patch.set_alpha(0.0)

    gs = GridSpec(1, len(labels) + 1,
                  width_ratios=[1] + [3] * len(labels), figure=fig)

    # Name labels column
    ax_names = fig.add_subplot(gs[0, 0])
    ax_names.patch.set_alpha(0.0)

    names = exon_info['name'].values
    unique_names = sorted(set(names))
    color_palette = plt.cm.tab10.colors[:len(unique_names)]
    name_colors = {n: color_palette[i] for i, n in enumerate(unique_names)}

    name_matrix = np.zeros((len(sorted_exon_ids), 1))
    name_cmap = LinearSegmentedColormap.from_list(
        'name_cmap', [(1, 1, 1)] + list(color_palette), N=len(unique_names) + 1
    )
    for i, name in enumerate(names):
        name_matrix[i, 0] = unique_names.index(name) + 1

    sns.heatmap(name_matrix, ax=ax_names, cmap=name_cmap, cbar=False,
                linewidths=0, rasterized=True)

    # Name group labels
    name_groups = {}
    current_name = None
    start_idx = 0
    for i, name in enumerate(names):
        if name != current_name:
            if current_name is not None:
                name_groups[current_name] = (start_idx, i - 1)
            current_name = name
            start_idx = i
    if current_name is not None:
        name_groups[current_name] = (start_idx, len(names) - 1)

    for name, (start, end) in name_groups.items():
        middle = (start + end) / 2
        ax_names.text(0.5, middle, name,
                      fontsize=10, fontweight='bold', ha='center', va='center',
                      color='black')

    ax_names.set_title('Name')
    ax_names.set_xticks([])
    ax_names.set_yticks([])

    # Plot each region
    for i, label in enumerate(labels):
        if label not in label_data:
            ax = fig.add_subplot(gs[0, i + 1])
            ax.set_facecolor('none')
            ax.text(0.5, 0.5, f"No data for {label}", ha='center', va='center')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(label)
            continue

        pivot = label_data[label]
        position_cols = sorted(pivot.columns)

        if '3ss' in label:
            min_pos, max_pos = 0, window + 50
        else:
            min_pos, max_pos = window - 50, window * 2

        position_cols = [pos for pos in position_cols if min_pos <= pos <= max_pos]

        display_matrix = np.zeros((len(sorted_exon_ids), len(position_cols)))
        for j, exon_id in enumerate(sorted_exon_ids):
            if exon_id in pivot.index:
                for k, pos in enumerate(position_cols):
                    if pos in pivot.columns:
                        display_matrix[j, k] = pivot.loc[exon_id, pos]

        ax = fig.add_subplot(gs[0, i + 1])
        ax.set_facecolor('none')
        sns.heatmap(display_matrix, ax=ax, cmap=colormaps['viridis'],
                    cbar=False, linewidths=0, rasterized=True)

        ax.set_title(label)
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.xaxis.set_visible(False)
        ax.set_yticks([])

        for name, (start, end) in name_groups.items():
            if end < len(sorted_exon_ids) - 1:
                ax.axhline(y=end + 1, color='white', linewidth=2, alpha=1)

    plt.tight_layout(rect=[0, 0, 0.95, 0.95])
    plt.savefig(f'{output_dir}/{FILEname}_heatmap.pdf',
                dpi=300, bbox_inches='tight')
    logging.info(f"Saved heatmap to {output_dir}/{FILEname}_heatmap.pdf")


def plot_rna_map(plotting_df, exon_categories, original_counts,
                 window, all_sites, output_dir, FILEname):
    """Generate the main RNA map -log10(pvalue) line plots."""
    sns.set(rc={'figure.figsize': (7, 5)})
    sns.set_style("whitegrid")

    if not all_sites:
        col_order = ["upstream_5ss", "middle_3ss", "middle_5ss", "downstream_3ss"]
        titles = ["Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS"]
        col_wrap = 4
    else:
        col_order = ["upstream_3ss", "upstream_5ss", "middle_3ss", "middle_5ss",
                     "downstream_3ss", "downstream_5ss"]
        titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS",
                  "Downstream 3'SS", "Downstream 5'SS"]
        col_wrap = 6

    g = sns.relplot(
        data=plotting_df, x='position', y='-log10pvalue_smoothed',
        hue='name', col='label', facet_kws={"sharex": False},
        kind='line', col_wrap=col_wrap, height=5, aspect=4 / 5,
        col_order=col_order
    )

    for ax, title in zip(g.axes.flat, titles):
        ax.set_title(title)
        ax.axhline(y=0, color='k', alpha=0.2, linewidth=0.5)
        fig = plt.gcf()
        marker_ax = add_enrichment_marker(fig, ax)

    g.set(xlabel='')
    g.axes[0].set_ylabel('-log10(p value) enrichment / control')

    sns.move_legend(
        g, "upper right",
        bbox_to_anchor=(1, 2),
        ncol=1, title=None, frameon=False
    )
    leg = g._legend
    set_legend_text(leg, exon_categories, original_counts)

    # Exon-intron drawings
    rect_fraction = 1 / ((window + 50) / 50)

    for i, ss_type in enumerate(col_order):
        ax = g.axes[i]
        is_middle = ss_type.startswith('middle_')
        exon_color = "midnightblue" if is_middle else "slategrey"

        if ss_type.endswith('_3ss'):
            ax.set_xlim([0, window + 50])
            ticks = np.arange(0, window + 51, 50)
            labels = ["" if t in (ticks[0], ticks[-1])
                      else str(int(t - window)) for t in ticks]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color=exon_color, alpha=1,
                transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)
        else:
            ax.set_xlim([window - 50, window * 2])
            ticks = np.arange(window - 50, window * 2 + 1, 50)
            labels = ["" if t in (ticks[0], ticks[-1])
                      else str(int(t - window)) for t in ticks]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color=exon_color, alpha=1,
                transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)

    plt.subplots_adjust(wspace=0.05)
    plt.savefig(
        f'{output_dir}/{FILEname}_RNAmap_-log10pvalue.pdf',
        bbox_extra_artists=([leg, rect, marker_ax]),
        bbox_inches='tight', pad_inches=0.8
    )
    logging.info(f"Saved RNA map to {output_dir}/{FILEname}_RNAmap_-log10pvalue.pdf")
    pbt.helpers.cleanup()


def plot_multivalency(middle_3ss_bed, middle_5ss_bed, downstream_3ss_bed,
                      upstream_5ss_bed, downstream_5ss_bed, upstream_3ss_bed,
                      fai, window, genome_fasta, output_dir, FILEname,
                      germsdir, all_sites, exon_categories, original_counts):
    """Run multivalency analysis and generate plots."""
    rect_fraction = 1 / ((window + 50) / 50)

    middle_3ss_mdf = get_multivalency_scores(
        middle_3ss_bed, fai, window, genome_fasta, output_dir,
        FILEname, 'middle_3ss', germsdir)
    middle_5ss_mdf = get_multivalency_scores(
        middle_5ss_bed, fai, window, genome_fasta, output_dir,
        FILEname, 'middle_5ss', germsdir)
    downstream_3ss_mdf = get_multivalency_scores(
        downstream_3ss_bed, fai, window, genome_fasta, output_dir,
        FILEname, 'downstream_3ss', germsdir)
    upstream_5ss_mdf = get_multivalency_scores(
        upstream_5ss_bed, fai, window, genome_fasta, output_dir,
        FILEname, 'upstream_5ss', germsdir)

    if all_sites:
        downstream_5ss_mdf = get_multivalency_scores(
            downstream_5ss_bed, fai, window, genome_fasta, output_dir,
            FILEname, 'downstream_5ss', germsdir)
        upstream_3ss_mdf = get_multivalency_scores(
            upstream_3ss_bed, fai, window, genome_fasta, output_dir,
            FILEname, 'upstream_3ss', germsdir)

    a = middle_3ss_mdf[0]
    b = middle_5ss_mdf[0]
    c = downstream_3ss_mdf[0]
    f = upstream_5ss_mdf[0]

    if not all_sites:
        plotting_df = pd.concat([f, a, b, c])
        mv_col_order = ["upstream_5ss", "middle_3ss", "middle_5ss", "downstream_3ss"]
        mv_titles = ["Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS"]
        mv_col_wrap = 4
    else:
        d = downstream_5ss_mdf[0]
        e = upstream_3ss_mdf[0]
        plotting_df = pd.concat([a, b, c, d, e, f])
        mv_col_order = ["upstream_3ss", "upstream_5ss", "middle_3ss", "middle_5ss",
                        "downstream_3ss", "downstream_5ss"]
        mv_titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS",
                     "Downstream 3'SS", "Downstream 5'SS"]
        mv_col_wrap = 6

    plotting_df.to_csv(
        f'{output_dir}/{FILEname}_RNAmap_multivalency.tsv', sep="\t"
    )
    logging.info(f"Saved multivalency TSV")

    # --- Main multivalency plot ---
    plt.figure()
    sns.set_style("whitegrid")

    g = sns.relplot(
        data=plotting_df, x='position', y='smoothed_kmer_multivalency',
        hue='exon_type', col='type', facet_kws={"sharex": False},
        kind='line', col_wrap=mv_col_wrap, height=5, aspect=3.5 / 5,
        errorbar=None, col_order=mv_col_order
    )

    for ax, title in zip(g.axes.flat, mv_titles):
        ax.set_title(title)
    g.set(xlabel='')
    g.axes[0].set_ylabel('mean smoothed kmer multivalency')
    leg = g._legend
    set_legend_text(leg, exon_categories, original_counts)

    _apply_mv_axis_formatting(g, mv_col_order, window, rect_fraction)

    plt.subplots_adjust(wspace=0.01)
    plt.savefig(
        f'{output_dir}/{FILEname}_RNAmap_multivalency.pdf',
        bbox_extra_artists=([leg, rect_fraction]),
        bbox_inches='tight', pad_inches=0.5
    )
    pbt.helpers.cleanup()
    logging.info(f"Saved multivalency plot")

    # --- Kmer plots ---
    a = middle_3ss_mdf[1]
    b = middle_5ss_mdf[1]
    c = downstream_3ss_mdf[1]
    f = upstream_5ss_mdf[1]

    if not all_sites:
        plotting_df = pd.concat([f, a, b, c])
    else:
        d = downstream_5ss_mdf[1]
        e = upstream_3ss_mdf[1]
        plotting_df = pd.concat([a, b, c, d, e, f])

    plotting_df.to_csv(
        f'{output_dir}/{FILEname}_RNAmap_TOP10KMER_multivalency.tsv', sep="\t"
    )

    # Silenced kmer plot
    _plot_kmer_multivalency(
        plotting_df[plotting_df['exon_type'] == 'silenced'],
        mv_col_order, mv_titles, mv_col_wrap, window, rect_fraction,
        output_dir, FILEname, 'silenced'
    )

    # Enhanced kmer plot
    _plot_kmer_multivalency(
        plotting_df[plotting_df['exon_type'] == 'enhanced'],
        mv_col_order, mv_titles, mv_col_wrap, window, rect_fraction,
        output_dir, FILEname, 'enhanced'
    )


def _apply_mv_axis_formatting(g, mv_col_order, window, rect_fraction):
    """Apply exon-intron axis formatting for multivalency plots."""
    for i, ss_type in enumerate(mv_col_order):
        ax = g.axes[i]
        if i == 0:
            ax.set_ylim(ymin=1)

        is_middle = ss_type.startswith('middle_')
        exon_color = "midnightblue" if is_middle else "slategrey"

        if ss_type.endswith('_3ss'):
            ax.set_xlim([window, (2 * window) + 50])
            ticks = np.arange(window, 2 * window + 50, 50)
            labels = ["" if t == ticks[0] else str(int(t - window)) for t in ticks]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color=exon_color, alpha=1, transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1, transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)
        else:
            ax.set_xlim([2 * window - 50, 3 * window])
            ticks = np.arange(2 * window - 50, 3 * window, 50)
            labels = ["" if t == ticks[0] else str(int(t - 2 * window)) for t in ticks]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color=exon_color, alpha=1, transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1, transform=ax.transAxes, clip_on=False)
            ax.add_artist(rect)


def _plot_kmer_multivalency(filtered_df, mv_col_order, mv_titles,
                            mv_col_wrap, window, rect_fraction,
                            output_dir, FILEname, exon_type):
    """Plot kmer multivalency for a given exon type."""
    plt.figure()
    sns.set_style("whitegrid")
    g = sns.relplot(
        data=filtered_df, x='position', y='smoothed_kmer_multivalency',
        hue='kmer', col='type', facet_kws={"sharex": False},
        kind='line', col_wrap=mv_col_wrap, height=5, errorbar=None,
        aspect=3.5 / 5, col_order=mv_col_order
    )
    for ax, title in zip(g.axes.flat, mv_titles):
        ax.set_title(title)
    g.set(xlabel='')
    g.axes[0].set_ylabel('mean smoothed kmer multivalency')
    leg = g._legend

    _apply_mv_axis_formatting(g, mv_col_order, window, rect_fraction)

    plt.subplots_adjust(wspace=0.01)
    plt.savefig(
        f'{output_dir}/{FILEname}_RNAmap_{exon_type}KMER_multivalency.pdf',
        bbox_extra_artists=([leg]),
        bbox_inches='tight', pad_inches=0.5
    )
    pbt.helpers.cleanup()
    logging.info(f"Saved {exon_type} kmer multivalency plot")


# ===============================================================================
# MAIN PIPELINE
# ===============================================================================

def run_rna_map(args):
    """
    Main RNA map pipeline. Handles both input modes with shared downstream logic.
    """
    output_dir = args.outputpath
    os.makedirs(output_dir, exist_ok=True)

    log_filename, start_time, logger = setup_logging(output_dir)
    logging.info(f"Log file created: {log_filename}")
    logging.info(f"Arguments: {args}")

    try:
        # Set random seed for reproducibility
        np.random.seed(args.seed)
        logging.info(f"Random seed set to {args.seed}")

        # Load chromosome list
        df_fai = pd.read_csv(args.fastaindex, sep='\t', header=None)
        chroms = set(df_fai[0].values)

        # ==============================================================
        # MODE SELECTION: Two tracks, one output format
        # ==============================================================

        if args.inputsplice:
            # ----- TRACK 1: rMATS -----
            input_mode = 'rmats'
            df_rmats = load_rmats_data(
                args.inputsplice,
                args.minctrl, args.maxctrl, args.maxincl,
                args.maxfdr, args.maxenh, args.minsil,
                chroms, args.no_constitutive
            )

            if args.prefix:
                FILEname = (args.prefix + "_" +
                            args.inputsplice.split('/')[-1]
                            .replace('.txt', '').replace('.gz', ''))
            else:
                FILEname = (args.inputsplice.split('/')[-1]
                            .replace('.txt', '').replace('.gz', ''))

        else:
            # ----- TRACK 2: VastDB -----
            input_mode = 'vastdb'
            df_rmats = load_vastdb_data(
                args.vastdb_enhanced,
                args.vastdb_silenced,
                args.vastdb_control,
                args.vastdb_constitutive,
                args.vastdb_annotation,
                chroms
            )

            if args.prefix:
                FILEname = args.prefix
            else:
                FILEname = "VastDB_RNAmap"

        # ==============================================================
        # SHARED PIPELINE: Same for both modes from here on
        # ==============================================================

        # Filter to valid chromosomes
        df_rmats = df_rmats[df_rmats['chr'].isin(chroms)]

        # Remove exons with missing flanking coordinates
        logging.info("\nFiltering exons with complete flanking coordinates...")
        before = len(df_rmats)
        df_rmats = df_rmats.dropna(
            subset=['upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
        )
        after = len(df_rmats)
        logging.info(f"Removed {before - after} exons with missing flanking coordinates")
        logging.info(f"Remaining: {after} exons")

        exon_categories = df_rmats.groupby('category').size()
        logging.info("\nExons in each category:")
        logging.info(exon_categories)

        # Validate categories
        if "control" not in exon_categories or exon_categories.loc["control"] == 0:
            logging.error("No control exons found!")
            sys.exit(1)

        if ("enhanced" not in exon_categories
                and "silenced" not in exon_categories):
            logging.error("No regulated exons found!")
            sys.exit(1)

        # Apply subsetting
        if not args.no_subset:
            df_rmats, original_counts = apply_subsetting(
                df_rmats, args.no_constitutive
            )
        else:
            logging.info("Subsetting disabled (--no_subset flag)")
            category_counts = df_rmats['category'].value_counts()
            original_counts = {cat: count for cat, count in category_counts.items()}

        exon_categories = df_rmats.groupby('category').size()

        # Save categorised exons
        if input_mode == 'rmats':
            save_cols = ['chr', 'exonStart_0base', 'exonEnd', 'strand', 'category',
                         'FDR', 'dPSI', 'maxPSI',
                         'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
            save_cols = [c for c in save_cols if c in df_rmats.columns]
            suffix = '_RMATS_with_categories.tsv'
        else:
            save_cols = ['chr', 'exonStart_0base', 'exonEnd', 'strand', 'category',
                         'EVENT', 'GENE',
                         'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
            save_cols = [c for c in save_cols if c in df_rmats.columns]
            suffix = '_VastDB_with_categories.tsv'

        df_rmats[save_cols].to_csv(
            f'{output_dir}/{FILEname}{suffix}', sep="\t", index=False
        )
        logging.info(f"Saved categorised exons to {output_dir}/{FILEname}{suffix}")

        # Exon length plots
        plot_exon_lengths(df_rmats.copy(), output_dir, FILEname)

        # ==============================================================
        # CREATE BED FILES FOR 6 SPLICE SITE REGIONS
        # ==============================================================
        # After strand correction at load time, upstream/downstream are
        # always in transcript order. Same-exon pairing is correct for
        # both modes.
        # ==============================================================
        logging.info("\n" + "=" * 60)
        logging.info("CREATING BED FILES FOR SPLICE SITES")
        logging.info("=" * 60)

        middle_3ss_bed = get_ss_bed(df_rmats, 'exonStart_0base', 'exonEnd')
        middle_5ss_bed = get_ss_bed(df_rmats, 'exonEnd', 'exonStart_0base')
        downstream_3ss_bed = get_ss_bed(df_rmats, 'downstreamES', 'downstreamEE')
        upstream_5ss_bed = get_ss_bed(df_rmats, 'upstreamEE', 'upstreamES')

        downstream_5ss_bed = None
        upstream_3ss_bed = None
        if args.all_sites:
            downstream_5ss_bed = get_ss_bed(df_rmats, 'downstreamEE', 'downstreamES')
            upstream_3ss_bed = get_ss_bed(df_rmats, 'upstreamES', 'upstreamEE')

        # ==============================================================
        # CALCULATE COVERAGE
        # ==============================================================
        if args.inputxlsites is not None:
            logging.info("\n" + "=" * 60)
            logging.info("CALCULATING COVERAGE")
            logging.info("=" * 60)

            fai = args.fastaindex
            xl_bed = args.inputxlsites
            window = args.window
            smoothing = args.smoothing

            middle_3ss = get_coverage_plot(
                xl_bed, middle_3ss_bed, fai, window, exon_categories,
                'middle_3ss', smoothing)
            middle_5ss = get_coverage_plot(
                xl_bed, middle_5ss_bed, fai, window, exon_categories,
                'middle_5ss', smoothing)
            downstream_3ss = get_coverage_plot(
                xl_bed, downstream_3ss_bed, fai, window, exon_categories,
                'downstream_3ss', smoothing)
            upstream_5ss = get_coverage_plot(
                xl_bed, upstream_5ss_bed, fai, window, exon_categories,
                'upstream_5ss', smoothing)

            linegraph_middle_3ss = middle_3ss[0]
            linegraph_middle_5ss = middle_5ss[0]
            linegraph_downstream_3ss = downstream_3ss[0]
            linegraph_upstream_5ss = upstream_5ss[0]

            heatmap_middle_3ss = middle_3ss[1]
            heatmap_middle_5ss = middle_5ss[1]
            heatmap_downstream_3ss = downstream_3ss[1]
            heatmap_upstream_5ss = upstream_5ss[1]

            if not args.all_sites:
                plotting_df = pd.concat([
                    linegraph_upstream_5ss, linegraph_middle_3ss,
                    linegraph_middle_5ss, linegraph_downstream_3ss
                ])
                heat_df = pd.concat([
                    heatmap_upstream_5ss, heatmap_middle_3ss,
                    heatmap_middle_5ss, heatmap_downstream_3ss
                ])
            else:
                downstream_5ss = get_coverage_plot(
                    xl_bed, downstream_5ss_bed, fai, window,
                    exon_categories, 'downstream_5ss', smoothing)
                upstream_3ss = get_coverage_plot(
                    xl_bed, upstream_3ss_bed, fai, window,
                    exon_categories, 'upstream_3ss', smoothing)

                linegraph_downstream_5ss = downstream_5ss[0]
                linegraph_upstream_3ss = upstream_3ss[0]
                heatmap_downstream_5ss = downstream_5ss[1]
                heatmap_upstream_3ss = upstream_3ss[1]

                plotting_df = pd.concat([
                    linegraph_middle_3ss, linegraph_middle_5ss,
                    linegraph_downstream_3ss, linegraph_downstream_5ss,
                    linegraph_upstream_3ss, linegraph_upstream_5ss
                ])
                heat_df = pd.concat([
                    heatmap_middle_3ss, heatmap_middle_5ss,
                    heatmap_downstream_3ss, heatmap_downstream_5ss,
                    heatmap_upstream_3ss, heatmap_upstream_5ss
                ])

            # Save coverage data
            plotting_df.to_csv(
                f'{output_dir}/{FILEname}_RNAmap.tsv', sep="\t", index=False
            )

            # Heatmap
            plot_heatmap(heat_df, exon_categories, window, args.all_sites,
                         output_dir, FILEname)

            # Main RNA map plot
            logging.info("\n" + "=" * 60)
            logging.info("PLOTTING RNA MAPS")
            logging.info("=" * 60)

            plot_rna_map(plotting_df, exon_categories, original_counts,
                         window, args.all_sites, output_dir, FILEname)

        # ==============================================================
        # MULTIVALENCY (optional, requires germs.R)
        # ==============================================================
        if hasattr(args, 'multivalency') and args.multivalency:
            logging.info("\n" + "=" * 60)
            logging.info("MULTIVALENCY ANALYSIS")
            logging.info("=" * 60)

            plot_multivalency(
                middle_3ss_bed, middle_5ss_bed,
                downstream_3ss_bed, upstream_5ss_bed,
                downstream_5ss_bed, upstream_3ss_bed,
                args.fastaindex, args.window, args.genomefasta,
                output_dir, FILEname, args.germsdir,
                args.all_sites, exon_categories, original_counts
            )

        logging.info("\n" + "=" * 60)
        logging.info("SCRIPT COMPLETED SUCCESSFULLY")
        logging.info("=" * 60)

    finally:
        log_runtime(start_time, logger)
        for handler in logger.handlers:
            handler.flush()
        logging.shutdown()


# ===============================================================================
# COMMAND LINE INTERFACE
# ===============================================================================

def cli():
    parser = argparse.ArgumentParser(
        description='Plot CLIP crosslinks around regulated exons to study '
                    'position-dependent impact on pre-mRNA splicing.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input modes:

  rMATS mode (original):
    python rna_maps.py -i rMATS.SE.MATS.JC.txt -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p PTBP1

  VastDB mode (ID lists):
    python rna_maps.py --vastdb_mode \\
      --vastdb_enhanced enhanced.txt --vastdb_silenced silenced.txt \\
      --vastdb_control control.txt --vastdb_constitutive constitutive.txt \\
      --vastdb_annotation EVENT_INFO-hg38.tab \\
      -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p AQR_K562
        """
    )

    # INPUT MODE - mutually exclusive
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-i', '--inputsplice', type=str,
        help='rMATS differential splicing file (rMATS mode)'
    )
    input_group.add_argument(
        '--vastdb_mode', action='store_true',
        help='Use VastDB ID lists mode (requires --vastdb_* arguments)'
    )

    # VASTDB-SPECIFIC ARGUMENTS
    vastdb_group = parser.add_argument_group('VastDB mode options')
    vastdb_group.add_argument(
        '--vastdb_enhanced', help='Enhanced exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_silenced', help='Silenced exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_control', help='Control exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_constitutive', help='Constitutive exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_annotation', help='VastDB EVENT_INFO file (e.g. EVENT_INFO-hg38.tab)')

    # SHARED REQUIRED ARGUMENTS
    required = parser.add_argument_group('Required arguments (both modes)')
    required.add_argument(
        '-x', '--inputxlsites', type=str, nargs='?',
        help='CLIP crosslinks in BED file format')
    required.add_argument(
        '-f', '--genomefasta', type=str, required=True,
        help='Genome FASTA file (.fa)')
    required.add_argument(
        '-fi', '--fastaindex', type=str, required=True,
        help='Genome FASTA index (.fai)')

    # SHARED OPTIONAL ARGUMENTS
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument(
        '-o', '--outputpath', type=str, default=os.getcwd(), nargs='?',
        help='Output folder [DEFAULT: current directory]')
    optional.add_argument(
        '-w', '--window', type=int, default=300, nargs='?',
        help='Window around splice sites [DEFAULT: 300]')
    optional.add_argument(
        '-s', '--smoothing', type=int, default=15, nargs='?',
        help='Smoothing window [DEFAULT: 15]')
    optional.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducible subsetting [DEFAULT: 42]')
    optional.add_argument(
        '-nc', '--no_constitutive', action='store_true',
        help='Exclude constitutive category')
    optional.add_argument(
        '-ns', '--no_subset', action='store_true',
        help='Disable subsetting of control/constitutive exons')
    optional.add_argument(
        '-ao', '--all_sites', action='store_true',
        help='Include all 6 splice sites (default: 4 core sites)')
    optional.add_argument(
        '-p', '--prefix', type=str,
        help='Prefix for output files')

    # rMATS-SPECIFIC THRESHOLDS
    rmats_group = parser.add_argument_group('rMATS mode thresholds')
    rmats_group.add_argument(
        '-mc', '--minctrl', type=float, default=-0.05, nargs='?',
        help='Minimum dPSI for control events [DEFAULT: -0.05]')
    rmats_group.add_argument(
        '-xc', '--maxctrl', type=float, default=0.05, nargs='?',
        help='Maximum dPSI for control events [DEFAULT: 0.05]')
    rmats_group.add_argument(
        '-xi', '--maxincl', type=float, default=0.9, nargs='?',
        help='Maximum PSI for control (above = constitutive) [DEFAULT: 0.9]')
    rmats_group.add_argument(
        '-xf', '--maxfdr', type=float, default=0.1, nargs='?',
        help='Maximum FDR for regulated events [DEFAULT: 0.1]')
    rmats_group.add_argument(
        '-xe', '--maxenh', type=float, default=-0.05, nargs='?',
        help='Maximum dPSI for enhanced exons [DEFAULT: -0.05]')
    rmats_group.add_argument(
        '-ms', '--minsil', type=float, default=0.05, nargs='?',
        help='Minimum dPSI for silenced exons [DEFAULT: 0.05]')

    # MULTIVALENCY (rMATS mode)
    mv_group = parser.add_argument_group('Multivalency analysis')
    mv_group.add_argument(
        '-v', '--multivalency', action='store_true',
        help='Run multivalency analysis (requires germs.R)')
    mv_group.add_argument(
        '-g', '--germsdir', type=str, default=os.getcwd(), nargs='?',
        help='Directory containing germs.R [DEFAULT: current directory]')

    args = parser.parse_args()

    # Validate VastDB mode requirements
    if args.vastdb_mode:
        if not args.vastdb_annotation:
            parser.error("--vastdb_mode requires --vastdb_annotation "
                         "(EVENT_INFO file)")
        if not any([args.vastdb_enhanced, args.vastdb_silenced,
                    args.vastdb_control, args.vastdb_constitutive]):
            parser.error("--vastdb_mode requires at least one "
                         "--vastdb_* ID list file")

    return args


# ===============================================================================
# ENTRY POINT
# ===============================================================================

if __name__ == '__main__':
    args = cli()
    run_rna_map(args)
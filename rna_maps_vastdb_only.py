#!/usr/bin/env python3
"""
RNA Maps from VastDB ID Lists

Generates RNA maps showing RBP binding patterns around categorized exons.
Uses VastDB EVENT IDs with pre-assigned categories (enhanced, silenced, control, constitutive).

Key features:
- 6 splice site regions (upstream_3ss, upstream_5ss, middle_3ss, middle_5ss, downstream_3ss, downstream_5ss)
- Coverage plots with statistical enrichment
- Heatmaps
- No rMATS dependency
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

sns.set_style("whitegrid", {'legend.frameon':True})
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
# VASTDB PARSING FUNCTIONS
# ===============================================================================

def parse_event_info_coordinates(event_info_df):
    """
    Parse genomic coordinates from EVENT_INFO file.
    Returns DataFrame with coordinates and flanking exon info.
    """
    coords_list = []
    
    for idx, row in event_info_df.iterrows():
        # Parse main exon coordinates from COORD_o
        coord_str = str(row['COORD_o'])
        match = re.match(r'(chr[\w]+):(\d+)-(\d+)', coord_str)
        if not match:
            continue
        
        # Extract strand from REF_CO
        ref_co = str(row.get('REF_CO', ''))
        strand_match = re.search(r':([+-])$', ref_co)
        strand = strand_match.group(1) if strand_match else '+'
        
        # Parse flanking exon coordinates from CO_C1 and CO_C2
        # CO_C1 format: "chr:start-end" (upstream constitutive exon)
        # CO_C2 format: "chr:start-end" (downstream constitutive exon)
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


def load_vastdb_categories(enhanced_file, silenced_file, control_file, 
                           constitutive_file, event_info_file):
    """
    Load VastDB ID lists and look up coordinates.
    Returns DataFrame ready for RNA map generation.
    """
    
    logging.info("="*60)
    logging.info("LOADING VASTDB ID LISTS")
    logging.info("="*60)
    
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
        id_to_category[eid] = 'constituitive'  # Match original script spelling
    
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
    logging.info(f"Matched {n_matched}/{n_total} IDs to coordinates ({100*n_matched/n_total:.1f}%)")
    
    if n_matched == 0:
        raise ValueError("No EVENT IDs matched to coordinates!")
    
    # Adjust coordinates from 1-based to 0-based
    event_coords_subset['exonStart_0base'] = event_coords_subset['exonStart_1based'] - 1
    
    # Also adjust flanking exon coordinates from 1-based to 0-based
    event_coords_subset['upstreamES'] = event_coords_subset['upstreamES'] - 1
    event_coords_subset['downstreamES'] = event_coords_subset['downstreamES'] - 1
    
    # Assign categories
    event_coords_subset['category'] = event_coords_subset['EVENT'].map(id_to_category)
    
    # Add fake FDR column (not used but needed for compatibility)
    event_coords_subset['FDR'] = 0.001
    
    # Create exon_id for tracking
    event_coords_subset['exon_id'] = (
        event_coords_subset['category'] + "_" + 
        event_coords_subset['chr'].astype(str) + ":" +
        event_coords_subset['exonStart_0base'].astype(str) + "-" +
        event_coords_subset['exonEnd'].astype(str) + ";" +
        event_coords_subset['strand'].astype(str)
    )
    
    logging.info("\nCategory distribution:")
    logging.info(event_coords_subset['category'].value_counts())
    
    return event_coords_subset


# ===============================================================================
# RNA MAP GENERATION FUNCTIONS (from original script)
# ===============================================================================

def get_ss_bed(df, pos_col, neg_col):
    """Create BED file for splice sites (handles strand orientation)."""
    # Strip "chr" prefix to match CLIP file format
    df = df.copy()

    ss_pos = df.loc[df['strand'] == "+", ['chr', pos_col, pos_col, 'exon_id', 'FDR', 'strand']]
    ss_pos.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_pos.start = ss_pos.start.transform(lambda x: x-1)

    ss_n = df.loc[df['strand'] == "-", ['chr', neg_col, neg_col, 'exon_id', 'FDR', 'strand']]
    ss_n.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_n.end = ss_n.end.transform(lambda x: x+1)

    ss = pd.concat([ss_pos, ss_n])

    return ss


def get_coverage_plot(xl_bed, df, fai, window, exon_categories, label, smoothing=15):
    """Calculate coverage of CLIP around splice sites with statistical enrichment."""
    df = df.loc[df.name != "."]
    xl_bed = pbt.BedTool(xl_bed).sort()
    pbt_df = pbt.BedTool.from_dataframe(
        df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    ).sort().slop(l=window, r=window, s=True, g=fai)
    
    df_coverage = pbt_df.coverage(
        b=xl_bed, **{'sorted': True, 's': True, 'd': True, 'nonamecheck': True}
    ).to_dataframe()[['thickStart', 'thickEnd', 'strand', 'name']]
    
    df_coverage.rename(columns=({'thickStart':'position','thickEnd':'coverage'}), inplace=True)

    df_plot = df_coverage
    
    # Adjust positions based on strand
    df_plot.loc[df_plot.strand=='+', 'position'] = df_plot['position'].astype('int32')
    df_plot.loc[df_plot.strand=='-', 'position'] = abs(2 * window + 2 - df_plot['position'])

    df_plot = df_plot.loc[df_plot.name != "."]

    df_plot = df_plot.join(
        df_plot.pop('name').str.split('_',n=1,expand=True).rename(columns={0:'name', 1:'exon_id'})
    )

    heatmap_plot = df_plot.copy()
    heatmap_plot['label'] = label

    # Aggregate coverage
    df_plot = df_plot.groupby(['name','position'], as_index=False).agg({'coverage':'sum'})
    
    exon_cat = pd.DataFrame({'name':exon_categories.index, 'number_exons':exon_categories.values})
    df_plot = df_plot.merge(exon_cat, how="left")

    df_plot['norm_coverage'] = np.where(
        df_plot['coverage'] == 0, 
        df_plot['coverage'], 
        df_plot['coverage']/df_plot['number_exons']
    )

    # Fisher test against control
    df_plot_ctrl = df_plot.loc[df_plot.name == "control"][["position","coverage","number_exons"]]
    df_plot_ctrl.columns = ["position","control_coverage","control_number_exons"]
    df_plot = df_plot.merge(df_plot_ctrl, how="left")
    df_plot['control_norm_coverage'] = df_plot["control_coverage"] / df_plot["control_number_exons"]
    
    # Pseudo-count for zero values
    df_plot.loc[df_plot['control_norm_coverage'] == 0, ['control_norm_coverage']] = 0.000001
    df_plot['fold_change'] = df_plot["norm_coverage"] / df_plot["control_norm_coverage"]
    
    # Statistical test
    contingency_table = list(zip(
        df_plot['coverage'], 
        df_plot['number_exons']-df_plot['coverage'], 
        df_plot['control_coverage'], 
        df_plot['control_number_exons'] - df_plot['control_coverage']
    ))
    contingency_table = [np.array(table).reshape(2,2) for table in contingency_table]
    df_plot['pvalue'] = [stats.fisher_exact(table)[1] for table in contingency_table]

    df_plot['-log10pvalue'] = np.log10(1/df_plot['pvalue'])
    df_plot['label'] = label

    df_plot.loc[df_plot['fold_change'] < 1, ['-log10pvalue']] = df_plot['-log10pvalue'] * -1
    df_plot['-log10pvalue_smoothed'] = df_plot['-log10pvalue'].rolling(
        smoothing, center=True, win_type="gaussian"
    ).mean(std=2)

    return df_plot, heatmap_plot


def set_legend_text(legend, exon_categories, original_counts=None):
    """Set legend text with optional subset information."""
    exon_cat = pd.DataFrame({'name': exon_categories.index, 'number_exons': exon_categories.values})
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
    y_zero_position = 0
    marker_height = 0.4
    normalized_bottom = (y_zero_position - y_min) / (y_max - y_min) - marker_height / 2
    
    marker_ax = fig.add_axes([0.94, normalized_bottom, 0.06, marker_height], 
                            label=f"enrichment_marker_{ax.get_title()}")
    
    marker_ax.set_xticks([])
    marker_ax.set_yticks([])
    for spine in marker_ax.spines.values():
        spine.set_visible(False)
    
    marker_ax.text(0.5, 0.5, "↑\nenriched\n\n\n↓\ndepleted",
                   ha='center', va='center', fontsize=11, 
                   transform=marker_ax.transAxes)
    
    return marker_ax


# ===============================================================================
# MAIN RNA MAP FUNCTION
# ===============================================================================

def run_rna_maps(enhanced_file, silenced_file, control_file, constitutive_file,
                event_info_file, xl_bed, fai, genome_fasta, window, smoothing,
                output_dir, no_constitutive, no_subset, all_sites, prefix):
    """
    Main function to generate RNA maps from VastDB ID lists.
    """
    
    log_filename, start_time, logger = setup_logging(output_dir)
    
    FILEname = prefix
    
    # Load chromosome list
    df_fai = pd.read_csv(fai, sep='\t', header=None)
    chroms = set(df_fai[0].values)
    
    # Load VastDB categories and coordinates
    df_rmats = load_vastdb_categories(
        enhanced_file, silenced_file, control_file, 
        constitutive_file, event_info_file
    )
    
    # Filter to valid chromosomes
    df_rmats = df_rmats[df_rmats['chr'].isin(chroms)]
    
    # Remove exons with missing flanking coordinates
    logging.info("\nFiltering exons with complete flanking coordinates...")
    before = len(df_rmats)
    df_rmats = df_rmats.dropna(subset=['upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'])
    after = len(df_rmats)
    logging.info(f"Removed {before - after} exons with missing flanking coordinates")
    logging.info(f"Remaining: {after} exons")
    
    exon_categories = df_rmats.groupby('category').size()
    logging.info("\nExons in each category:")
    logging.info(exon_categories)
    
    # Check for required categories
    if "control" not in exon_categories or exon_categories.loc["control"] == 0:
        raise ValueError("No control exons found!")
    
    if "enhanced" not in exon_categories and "silenced" not in exon_categories:
        raise ValueError("No regulated exons found!")
    
    # Apply subsetting
    if not no_subset:
        category_counts = df_rmats['category'].value_counts()
        original_counts = {cat: count for cat, count in category_counts.items()}

        np.random.seed(42) 
        logging.info("Random seed set to 42 for reproducible subsetting")
        
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
            control_indices_to_keep = np.random.choice(control_indices, target_count, replace=False)
            drop_mask = df_rmats.index.isin(control_indices) & ~df_rmats.index.isin(control_indices_to_keep)
            df_rmats = df_rmats[~drop_mask]
            logging.info(f"Subsetted control exons from {category_counts['control']} to {target_count}")
        
        # Subset constitutive
        if not no_constitutive and 'constituitive' in category_counts and category_counts['constituitive'] > target_count > 0:
            const_indices = df_rmats[df_rmats['category'] == 'constituitive'].index
            const_indices_to_keep = np.random.choice(const_indices, target_count, replace=False)
            drop_mask = df_rmats.index.isin(const_indices) & ~df_rmats.index.isin(const_indices_to_keep)
            df_rmats = df_rmats[~drop_mask]
            logging.info(f"Subsetted constitutive exons from {category_counts['constituitive']} to {target_count}")
    else:
        logging.info("Subsetting disabled")
        category_counts = df_rmats['category'].value_counts()
        original_counts = {cat: count for cat, count in category_counts.items()}
    
    exon_categories = df_rmats.groupby('category').size()
    
    # Save categorized exons
    output_file = os.path.join(output_dir, f"{FILEname}_VastDB_with_categories.tsv")
    df_rmats[['chr', 'exonStart_0base', 'exonEnd', 'strand', 'category', 'EVENT', 'GENE',
             'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']].to_csv(
        output_file, sep='\t', index=False
    )
    logging.info(f"Saved categorized exons to {output_file}")
    
    # Create BED files for 6 regions
    logging.info("\n" + "="*60)
    logging.info("CREATING BED FILES FOR SPLICE SITES")
    logging.info("="*60)
    
    middle_3ss_bed = get_ss_bed(df_rmats, 'exonStart_0base', 'exonEnd')
    middle_5ss_bed = get_ss_bed(df_rmats, 'exonEnd', 'exonStart_0base')
    downstream_3ss_bed = get_ss_bed(df_rmats, 'downstreamES', 'downstreamEE')
    upstream_5ss_bed = get_ss_bed(df_rmats, 'upstreamEE', 'upstreamES')
    
    if all_sites:
        downstream_5ss_bed = get_ss_bed(df_rmats, 'downstreamEE', 'downstreamES')
        upstream_3ss_bed = get_ss_bed(df_rmats, 'upstreamES', 'upstreamEE')
    
    # Calculate coverage
    if xl_bed is not None:
        logging.info("\n" + "="*60)
        logging.info("CALCULATING COVERAGE")
        logging.info("="*60)
        
        middle_3ss = get_coverage_plot(xl_bed, middle_3ss_bed, fai, window, exon_categories, 'middle_3ss', smoothing)
        middle_5ss = get_coverage_plot(xl_bed, middle_5ss_bed, fai, window, exon_categories, 'middle_5ss', smoothing)
        downstream_3ss = get_coverage_plot(xl_bed, downstream_3ss_bed, fai, window, exon_categories, 'downstream_3ss', smoothing)
        upstream_5ss = get_coverage_plot(xl_bed, upstream_5ss_bed, fai, window, exon_categories, 'upstream_5ss', smoothing)
        
        linegraph_middle_3ss = middle_3ss[0]
        linegraph_middle_5ss = middle_5ss[0]
        linegraph_downstream_3ss = downstream_3ss[0]
        linegraph_upstream_5ss = upstream_5ss[0]
        
        heatmap_middle_3ss = middle_3ss[1]
        heatmap_middle_5ss = middle_5ss[1]
        heatmap_downstream_3ss = downstream_3ss[1]
        heatmap_upstream_5ss = upstream_5ss[1]
        
        if not all_sites:
            plotting_df = pd.concat([linegraph_upstream_5ss, linegraph_middle_3ss, 
                                    linegraph_middle_5ss, linegraph_downstream_3ss])
            heat_df = pd.concat([heatmap_upstream_5ss, heatmap_middle_3ss, 
                                heatmap_middle_5ss, heatmap_downstream_3ss])
        else:
            downstream_5ss = get_coverage_plot(xl_bed, downstream_5ss_bed, fai, window, exon_categories, 'downstream_5ss', smoothing)
            upstream_3ss = get_coverage_plot(xl_bed, upstream_3ss_bed, fai, window, exon_categories, 'upstream_3ss', smoothing)
            
            linegraph_downstream_5ss = downstream_5ss[0]
            linegraph_upstream_3ss = upstream_3ss[0]
            heatmap_downstream_5ss = downstream_5ss[1]
            heatmap_upstream_3ss = upstream_3ss[1]
            
            plotting_df = pd.concat([linegraph_middle_3ss, linegraph_middle_5ss, 
                                    linegraph_downstream_3ss, linegraph_downstream_5ss, 
                                    linegraph_upstream_3ss, linegraph_upstream_5ss])
            heat_df = pd.concat([heatmap_middle_3ss, heatmap_middle_5ss, 
                                heatmap_downstream_3ss, heatmap_downstream_5ss, 
                                heatmap_upstream_3ss, heatmap_upstream_3ss])
        
        # Save data
        plotting_df.to_csv(f'{output_dir}/{FILEname}_RNAmap.tsv', sep="\t", index=False)
        
        # Plot RNA maps
        logging.info("\n" + "="*60)
        logging.info("PLOTTING RNA MAPS")
        logging.info("="*60)
        
        sns.set(rc={'figure.figsize':(7, 5)})
        sns.set_style("whitegrid")
        
        if not all_sites:
            col_order = ["upstream_5ss", "middle_3ss", "middle_5ss", "downstream_3ss"]
            titles = ["Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS"]
            col_wrap = 4
        else:
            col_order = ["upstream_3ss", "upstream_5ss", "middle_3ss", "middle_5ss", "downstream_3ss", "downstream_5ss"]
            titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS", "Downstream 5'SS"]
            col_wrap = 6
        
        g = sns.relplot(
            data=plotting_df, 
            x='position', 
            y='-log10pvalue_smoothed', 
            hue='name', 
            col='label', 
            facet_kws={"sharex":False},
            kind='line', 
            col_wrap=col_wrap, 
            height=5, 
            aspect=4/5,
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
            ncol=1,
            title=None,
            frameon=False
        )
        leg = g._legend
        set_legend_text(leg, exon_categories, original_counts)
        
        # Add exon-intron drawings
        rect_fraction = 1 / ((window + 50) / 50)
        
        for i, ss_type in enumerate(col_order):
            ax = g.axes[i]
            is_middle = ss_type.startswith('middle_')
            exon_color = "midnightblue" if is_middle else "slategrey"
            
            if ss_type.endswith('_3ss'):
                # 3' splice site: exon on right
                ax.set_xlim([0, window + 50])
                ticks = np.arange(0, window + 51, 50)
                labels = ["" if t in (ticks[0], ticks[-1]) else str(int(t - window)) for t in ticks]
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
                # 5' splice site: exon on left
                ax.set_xlim([window - 50, window * 2])
                ticks = np.arange(window - 50, window * 2 + 1, 50)
                labels = ["" if t in (ticks[0], ticks[-1]) else str(int(t - window)) for t in ticks]
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
        plt.savefig(f'{output_dir}/{FILEname}_RNAmap_-log10pvalue.pdf', 
                   bbox_extra_artists=([leg, rect, marker_ax]),
                   bbox_inches='tight',
                   pad_inches=0.8)
        
        logging.info(f"Saved RNA map to {output_dir}/{FILEname}_RNAmap_-log10pvalue.pdf")
        
        pbt.helpers.cleanup()
    
    log_runtime(start_time, logger)
    logging.info("="*60)
    logging.info("SCRIPT COMPLETED SUCCESSFULLY")
    logging.info("="*60)


# ===============================================================================
# COMMAND LINE INTERFACE
# ===============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Generate RNA maps from VastDB ID lists',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python rna_maps_vastdb_only.py \\
    --vastdb_enhanced enhanced_ids.txt \\
    --vastdb_silenced silenced_ids.txt \\
    --vastdb_control control_ids.txt \\
    --vastdb_constitutive constitutive_ids.txt \\
    --vastdb_annotation EVENT_INFO-hg38.tab \\
    -x CLIP.bed -f genome.fa -fi genome.fa.fai \\
    -o output -p PTBP1
        """
    )
    
    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--vastdb_annotation', required=True,
                         help='VastDB EVENT_INFO file')
    required.add_argument('-x', '--clip', required=True,
                         help='CLIP crosslinks BED file')
    required.add_argument('-f', '--fasta', required=True,
                         help='Genome FASTA file')
    required.add_argument('-fi', '--fasta_index', required=True,
                         help='Genome FASTA index (.fai)')
    required.add_argument('-o', '--output_dir', required=True,
                         help='Output directory')
    required.add_argument('-p', '--prefix', required=True,
                         help='Output file prefix')
    
    # VastDB ID lists
    vastdb_group = parser.add_argument_group('VastDB ID lists (at least one required)')
    vastdb_group.add_argument('--vastdb_enhanced',
                             help='Enhanced exon IDs (one per line)')
    vastdb_group.add_argument('--vastdb_silenced',
                             help='Silenced exon IDs (one per line)')
    vastdb_group.add_argument('--vastdb_control',
                             help='Control exon IDs (one per line)')
    vastdb_group.add_argument('--vastdb_constitutive',
                             help='Constitutive exon IDs (one per line)')
    
    # Optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-w', '--window', type=int, default=300,
                         help='Window around splice sites (default: 300)')
    optional.add_argument('-s', '--smoothing', type=int, default=15,
                         help='Smoothing window (default: 15)')
    optional.add_argument('-nc', '--no_constitutive', action='store_true',
                         help='Exclude constitutive exons')
    optional.add_argument('-ns', '--no_subset', action='store_true',
                         help='Disable subsetting')
    optional.add_argument('-ao', '--all_sites', action='store_true',
                         help='Include all 6 splice sites (default: 4 core sites)')
    
    args = parser.parse_args()
    
    # Validate at least one ID list provided
    if not any([args.vastdb_enhanced, args.vastdb_silenced, 
                args.vastdb_control, args.vastdb_constitutive]):
        parser.error("At least one VastDB ID list must be provided")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run RNA maps
    run_rna_maps(
        enhanced_file=args.vastdb_enhanced,
        silenced_file=args.vastdb_silenced,
        control_file=args.vastdb_control,
        constitutive_file=args.vastdb_constitutive,
        event_info_file=args.vastdb_annotation,
        xl_bed=args.clip,
        fai=args.fasta_index,
        genome_fasta=args.fasta,
        window=args.window,
        smoothing=args.smoothing,
        output_dir=args.output_dir,
        no_constitutive=args.no_constitutive,
        no_subset=args.no_subset,
        all_sites=args.all_sites,
        prefix=args.prefix
    )


if __name__ == '__main__':
    main()

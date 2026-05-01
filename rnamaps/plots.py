"""Main RNA map plots: exon length distributions, per-exon heatmap, RNA map line plot."""

import logging

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools as pbt
import seaborn as sns
from matplotlib import colormaps
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

from rnamaps.config import colors_dict
from rnamaps.plot_helpers import add_enrichment_marker, set_legend_text
from rnamaps.preprocessing import smooth_coverage


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
        order=["control", "constitutive", "enhanced", "enhanced_rest",
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
        row_index = pd.Index(sorted_exon_ids)
        col_index = pd.Index(position_cols)
        rows_present = row_index.intersection(pivot.index, sort=False)
        cols_present = col_index.intersection(pivot.columns, sort=False)

        if len(rows_present) and len(cols_present):
            row_pos = row_index.get_indexer(rows_present)
            col_pos = col_index.get_indexer(cols_present)
            display_matrix[np.ix_(row_pos, col_pos)] = pivot.loc[
                rows_present, cols_present
            ].to_numpy()

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

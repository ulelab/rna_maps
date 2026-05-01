"""Optional multivalency analysis using the Ule lab GeRMs R package."""

import logging
import os
import random
import string

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools as pbt
import seaborn as sns

from rnamaps.plot_helpers import set_legend_text


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

import pandas as pd
import numpy as np
import pybedtools as pbt
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import os
import random
import string
import logging

smoothing = 15

def categorize_exons(
    df,
    min_ctrl = -0.05, max_ctrl = 0.05,
    max_inclusion = 0.9,
    max_fdr = 0.1, max_enh = -0.05, min_sil = 0.05,
    no_constitutive = False
):

    df_cat = df.copy()

    # Define boolean masks
    sil_mask = (df_cat['dPSI'] > min_sil) & (df_cat['FDR'] < max_fdr)
    enh_mask = (df_cat['dPSI'] < max_enh) & (df_cat['FDR'] < max_fdr)
    ctrl_mask = (df_cat['dPSI'] > min_ctrl) & (df_cat['dPSI'] < max_ctrl)

    conditions = [sil_mask, enh_mask]
    choices = ['silenced', 'enhanced']

    if not no_constitutive:
        const_mask = ctrl_mask & (df_cat['maxPSI'] > max_inclusion)
        conditions.append(const_mask)
        choices.append('constitutive')
    # Control must be appended as last choice (even if constitutive exists)
    conditions.append(ctrl_mask)
    choices.append('control')

    df_cat['category'] = np.select(conditions, choices, default=None)
    return df_cat

def downsample_categories(df, seed = None):
    if seed is not None:
        np.random.seed(seed)

    df_down = df.copy()
    category_counts = df_down['category'].value_counts().to_dict()
    counts_before = category_counts.copy()

    # Determine target: max count among 'silenced' and 'enhanced'
    target = 0
    if 'enhanced' in category_counts and 'silenced' in category_counts:
        target = max(category_counts['enhanced'], category_counts['silenced'])
    elif 'enhanced' in category_counts:
        target = category_counts['enhanced']
    elif 'silenced' in category_counts:
        target = category_counts['silenced']

    # Helper to subsample a given category to target
    def subsample(cat_name):
        nonlocal df_down
        idx = df_down.index[df_down['category'] == cat_name].tolist()
        if len(idx) > target > 0:
            keep = np.random.choice(idx, size=target, replace=False)
            drop = set(idx) - set(keep)
            df_down = df_down.drop(index=list(drop))

    # Subsample 'control'
    if 'control' in category_counts:
        subsample('control')
    # Subsample 'constitutive' if exists
    if 'constitutive' in category_counts:
        subsample('constitutive')

    # Recompute counts_after
    counts_after = df_down['category'].value_counts().to_dict()
    return df_down.reset_index(drop=True), counts_after, counts_before

def compute_exon_lengths(df):
    df_len = df.copy()
    df_len['regulated_exon_length'] = df_len['exonEnd'] - df_len['exonStart_0base']
    df_len['upstream_exon_length'] = df_len['upstreamEE'] - df_len['upstreamES']
    df_len['downstream_exon_length'] = df_len['downstreamEE'] - df_len['downstreamES']

    length_df = df_len[[
        'category',
        'regulated_exon_length',
        'upstream_exon_length',
        'downstream_exon_length'
    ]]

    length_melted = length_df.melt(
        id_vars=['category'],
        var_name='exon_type',
        value_name='exon_length'
    )
    return length_melted

def smooth_coverage(df, window_size=10, std=2):
    # Create a copy of the dataframe to avoid modifying the original
    result = df.copy()

    # Process each unique exon_id and label combination
    groups = []
    for (exon_id, label), group_df in result.groupby(['exon_id', 'label']):
        # Skip if too few points to smooth
        if len(group_df) < 3:
            groups.append(group_df)
            continue

        # Sort by position for rolling window
        group_sorted = group_df.sort_values('position')

        # Apply smoothing only if we have enough points
        if len(group_sorted) >= window_size:
            # Get coverage values
            values = group_sorted['coverage'].values

            # Create a Series with integer indices (not the DataFrame indices)
            # This ensures clean math without index alignment issues
            s = pd.Series(values)

            # Apply rolling window
            smoothed = s.rolling(
                window=window_size,
                center=True,
                win_type='gaussian'
            ).mean(std=std)

            # Fill NaN values at the edges with original values
            smoothed = smoothed.fillna(s)

            # Update the coverage in the original group
            # This works because the sorted indices correspond to the smoothed series
            group_sorted['coverage'] = smoothed.values

        # Add the processed group to our list
        groups.append(group_sorted)

    # Combine all groups back into a single DataFrame
    return pd.concat(groups)

def get_coverage_plot(xl_bed, df, fai, window, exon_categories, label):
    """Return coverage of xl_bed items around df features extended by windows"""
    df = df.loc[df.name != "."]
    xl_bed = pbt.BedTool(xl_bed).sort()
    pbt_df = pbt.BedTool.from_dataframe(df[['chr', 'start', 'end', 'name', 'score', 'strand']]).sort().slop(l=window,
                                                                                                            r=window,
                                                                                                            s=True,
                                                                                                            g=fai)
    df_coverage = \
    pbt_df.coverage(b=xl_bed, **{'sorted': True, 's': True, 'd': True, 'nonamecheck': True}).to_dataframe()[
        ['thickStart', 'thickEnd', 'strand', 'name']]
    df_coverage.rename(columns=({'thickStart': 'position', 'thickEnd': 'coverage'}), inplace=True)

    df_plot = df_coverage

    df_plot.loc[df_plot.strand == '+', 'position'] = df_plot['position'].astype('int32')
    df_plot.loc[df_plot.strand == '-', 'position'] = abs(2 * window + 2 - df_plot['position'])

    df_plot = df_plot.loc[df_plot.name != "."]

    df_plot = df_plot.join(
        df_plot.pop('name').str.split('_', n=1, expand=True).rename(columns={0: 'name', 1: 'exon_id'}))

    heatmap_plot = df_plot
    heatmap_plot['label'] = label

    df_plot = df_plot.groupby(['name', 'position'], as_index=False).agg({'coverage': 'sum'})

    exon_cat = pd.DataFrame({'name': exon_categories.index, 'number_exons': exon_categories.values})
    df_plot = df_plot.merge(exon_cat, how="left")

    df_plot['norm_coverage'] = np.where(df_plot['coverage'] == 0, df_plot['coverage'],
                                        df_plot['coverage'] / df_plot['number_exons'])

    # fisher test [covered exons, non-covered exon, control covered exons, control non-covered exons]
    df_plot_ctrl = df_plot.loc[df_plot.name == "control"][["position", "coverage", "number_exons"]]
    df_plot_ctrl.columns = ["position", "control_coverage", "control_number_exons"]
    df_plot = df_plot.merge(df_plot_ctrl, how="left")
    df_plot['control_norm_coverage'] = df_plot["control_coverage"] / df_plot["control_number_exons"]
    # give 0 values a small pseudo-count so that fold change can be calculated
    df_plot.loc[df_plot['control_norm_coverage'] == 0, ['control_norm_coverage']] = 0.000001
    df_plot['fold_change'] = df_plot["norm_coverage"] / df_plot["control_norm_coverage"]

    contingency_table = list(
        zip(df_plot['coverage'], df_plot['number_exons'] - df_plot['coverage'], df_plot['control_coverage'],
            df_plot['control_number_exons'] - df_plot['control_coverage']))
    contingency_table = [np.array(table).reshape(2, 2) for table in contingency_table]
    df_plot['pvalue'] = [stats.fisher_exact(table)[1] for table in contingency_table]

    df_plot['-log10pvalue'] = np.log10(1 / df_plot['pvalue'])
    df_plot['label'] = label

    df_plot.loc[df_plot['fold_change'] < 1, ['-log10pvalue']] = df_plot['-log10pvalue'] * -1
    df_plot['-log10pvalue_smoothed'] = df_plot['-log10pvalue'].rolling(smoothing, center=True,
                                                                       win_type="gaussian").mean(std=2)

    return df_plot, heatmap_plot

def get_multivalency_scores(df, fai, window, genome_fasta, output_dir, name, type, germsdir):
    """Return multivalency scores around df features extended by windows"""
    df = df.loc[(df.name != ".") & (pd.notnull(df.name)) & (df.name != "None")]
    df = df[['chr', 'start', 'end', 'name', 'score', 'strand']]
    df.columns = ['chr', 'start', 'stop', 'name', 'score', 'strand']
    df['name'] = df['name'].apply(lambda x: (str(x) + "XX" + random.choice(list(string.ascii_lowercase)) + random.choice(list(string.ascii_lowercase)) + random.choice(list(string.ascii_lowercase)) + random.choice(list(string.ascii_lowercase))  + random.choice(list(string.ascii_lowercase)) + random.choice(list(string.ascii_lowercase))))
    pbt_df = pbt.BedTool.from_dataframe(df[['chr', 'start', 'stop', 'name', 'score', 'strand']]).sort().slop(l=2*window, r=2*window, s=True, g=fai)
    logging.info("Number of sites: " + str(len(pbt_df)))
    pbts = pbt.BedTool.filter(pbt_df, lambda x: len(x) == (4*window) + 1).saveas()
    logging.info("Number of seqs considered after filtering those that run off the end of chroms: " + str(len(pbts)))
    pbts.sequence(fi=genome_fasta,name=True).save_seqs(f'{output_dir}/{name}_{type}_temp.fa')
    logging.info("Running germs to calculate multivalency scores...")

    os.system("RScript --vanilla " + germsdir + "/germs.R -f " + f'{output_dir}/{name}_{type}_temp.fa' + " -w 100 -s 20")
    os.system("gunzip -f *multivalency.tsv.gz")
    mdf = pd.read_csv(f'{output_dir}/{name}_{type}_temp_5_101_21.multivalency.tsv', sep='\t', header=0)
    os.system(f'rm {output_dir}/{name}_{type}_temp_5_101_21.multivalency.tsv')
    os.system(f'rm {output_dir}/{name}_{type}_temp.fa')
    mdf['position'] = np.tile(np.arange(0, 4*window - 3), len(pbts))

    mdf[['exon_type','label','roname']] = mdf['sequence_name'].str.split(r'XX|_',expand=True)

    # GET MOST CONTRIBUTING KMERS
    # Filter for 'exon_type' being 'enhanced' or 'silenced'
    filtered_mdf = mdf[mdf['exon_type'].isin(['enhanced', 'silenced'])]
    # Separate dataframes for enhanced and silenced exon_types
    enhanced_df = filtered_mdf[filtered_mdf['exon_type'] == 'enhanced']
    silenced_df = filtered_mdf[filtered_mdf['exon_type'] == 'silenced']
    # Function to calculate top 10 kmers for a given DataFrame
    def calculate_top_kmers(df):
        kmer_scores = df.groupby('kmer')['smoothed_kmer_multivalency'].sum()
        top_kmers = kmer_scores.nlargest(5).index
        return df[df['kmer'].isin(top_kmers)]
    # Calculate top 10 kmers for enhanced and silenced exon_types separately
    top_kmers_enhanced = calculate_top_kmers(enhanced_df)
    top_kmers_silenced = calculate_top_kmers(silenced_df)

    # Combine the top kmers for enhanced and silenced into one table
    top_kmers_df = pd.concat([top_kmers_enhanced, top_kmers_silenced])

    mdf = mdf.groupby(['exon_type','position'], as_index=False).agg({'smoothed_kmer_multivalency':'mean'}).reset_index()
    top_kmers_df = top_kmers_df.groupby(['position','exon_type','kmer'], as_index=False).agg({'smoothed_kmer_multivalency':'mean'}).reset_index()

    mdf['type'] = type
    top_kmers_df['type'] = type

    print(mdf.head())
    print(top_kmers_df.head())

    # mdf = mdf.loc[(mdf.label != ".") & (pd.notnull(mdf.label)) & (mdf.label != "None")]
    # top_kmers_df = top_kmers_df.loc[(top_kmers_df.label != ".") & (pd.notnull(top_kmers_df.label)) & (top_kmers_df.label != "None")]

    return mdf,top_kmers_df

def do_multivalency_analysis(ss_map,fai,window,genome_fasta,output_dir,FILEname,germsdir):

    aggregated_list = []
    top10_list = []
    for label, bed_df in ss_map.items():
        agg_df, top10_df = get_multivalency_scores(
            bed_df, fai, window, genome_fasta, output_dir, FILEname, label, germsdir
        )
        # Ensure 'type' column in both DataFrames indicates label
        agg_df['type'] = label
        top10_df['type'] = label
        aggregated_list.append(agg_df)
        top10_list.append(top10_df)

    # Combine & save TSVs
    combined_agg = pd.concat(aggregated_list, ignore_index=True)
    combined_agg.to_csv(f"{output_dir}/{FILEname}_RNAmap_multivalency.tsv", sep="\t", index=False)

    combined_top10 = pd.concat(top10_list, ignore_index=True)
    combined_top10.to_csv(f"{output_dir}/{FILEname}_RNAmap_TOP10KMER_multivalency.tsv", sep="\t", index=False)

    # Plot aggregate multivalency curves
    plt.figure()
    sns.set_style("whitegrid")
    g1 = sns.relplot(
        data=combined_agg,
        x='position',
        y='smoothed_kmer_multivalency',
        hue='exon_type',
        col='type',
        kind='line',
        col_wrap=len(aggregated_list),
        height=4,
        aspect=1.2,
        facet_kws={"sharex": False},
        col_order=list(ss_map.keys()),
        errorbar=None
    )
    # (Customize subplots similarly to plot_enrichment_curves, adding exon schematic, etc.)
    plt.savefig(f"{output_dir}/{FILEname}_RNAmap_multivalency.pdf", dpi=300, bbox_inches='tight')
    plt.close()

    # Plot top10 k-mer curves for 'silenced' and 'enhanced'
    for category in ['silenced', 'enhanced']:
        subset = combined_top10[combined_top10['exon_type'] == category].copy()
        if subset.empty:
            continue
        plt.figure()
        sns.set_style("whitegrid")
        g2 = sns.relplot(
            data=subset,
            x='position',
            y='smoothed_kmer_multivalency',
            hue='kmer',
            col='type',
            kind='line',
            col_wrap=len(ss_map),
            height=4,
            aspect=1.2,
            facet_kws={"sharex": False},
            col_order=list(ss_map.keys()),
            errorbar=None
        )
        # (Customize subplots similarly)
        plt.savefig(f"{output_dir}/{FILEname}_RNAmap_{category}KMER_multivalency.pdf", dpi=300, bbox_inches='tight')
        plt.close()
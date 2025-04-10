import matplotlib
import matplotlib.ticker as mticker
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

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

def setup_logging(output_path):
    """Sets up logging to file and console, and returns the log filename + start time."""
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"execution_{timestamp}.log"
    log_path = os.path.join(output_path, log_filename)

    # Get root logger and reset handlers to avoid duplicates
    logger = logging.getLogger()
    logger.handlers.clear()
    logger.setLevel(logging.INFO)

    # Create file handler
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    # Create stream handler (console output)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    # Ensure file writes are flushed automatically
    file_handler.flush = file_handler.stream.flush

    # Capture start time
    start_time = time.time()

    logger.info(f"Starting script execution at {timestamp}")

    return log_filename, start_time, logger

def log_runtime(start_time, logger):
    """Logs the script execution time in minutes and seconds format."""
    end_time = time.time()
    total_seconds = int(end_time - start_time)

    minutes, seconds = divmod(total_seconds, 60)
    runtime_str = f"{minutes}m {seconds}s" if minutes > 0 else f"{seconds}s"

    # Ensure log is written immediately
    logger.info(f"Script completed in {runtime_str}.")
    for handler in logger.handlers:
        handler.flush()

sns.set_style("whitegrid", {'legend.frameon':True})

colors_dict = {'all': '#D3D3D3', 'ctrl': '#408F76', 'enh': '#F30C08', 'sil': '#005299', 'enhrest': '#FFB122', 'silrest': '#6DC2F5', 'const': '#666666'}
linewidth = 3
dashes = False

def cli():
    parser = argparse.ArgumentParser(description='Plot CLIP crosslinks around regulated exons to study position-dependent impact on pre-mRNA splicing.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i',"--inputsplice", type=str, required=True,
                        help='quantification of differential splicing produced by rMATS')
    optional.add_argument('-x',"--inputxlsites", type=str, nargs='?',
                        help='CLIP crosslinks in BED file format')
    required.add_argument('-f',"--genomefasta", type=str, required=True,
                        help='genome fasta file (.fa)')
    required.add_argument('-fi',"--fastaindex", type=str, required=True,
                        help='genome fasta index file (.fai)')
    optional.add_argument('-o',"--outputpath", type=str, default=os.getcwd(), nargs='?',
                        help='output folder [DEFAULT current directory]')
    optional.add_argument('-w',"--window", type=int, default=300, nargs='?',
                        help='window around regulated splicing events to plot crosslinks [DEFAULT 300]')
    optional.add_argument('-s',"--smoothing", type=int, default=15, nargs='?',
                        help='smoothing window for plotting crosslink signal [DEFAULT 15]')
    optional.add_argument('-mc',"--minctrl", type=float, default=-0.05, nargs='?',
                        help='minimum dPSI for control events [DEFAULT -0.05]')
    optional.add_argument('-xc',"--maxctrl", type=float, default=0.05, nargs='?',
                        help='maximum dPSI for control events [DEFAULT 0.05]')
    optional.add_argument('-xi',"--maxincl", type=float, default=0.9, nargs='?',
                        help='maximum PSI for control exons, above this limit exons are considered constitutive [DEFAULT 0.9]')
    optional.add_argument('-xf',"--maxfdr", type=float, default=0.1, nargs='?',
                        help='maximum FDR for regulated events, above this events fall in "rest" class, is used for rMATS [DEFAULT 0.1]')
    optional.add_argument('-xe',"--maxenh", type=float, default=-0.05, nargs='?',
                        help='maximum inclusion for exons to be considered enhanced [DEFAULT -0.05]')
    optional.add_argument('-ms',"--minsil", type=float, default=0.05, nargs='?',
                        help='minimum inclusion for exons to be considered silenced [DEFAULT 0.05]')
    optional.add_argument('-v','--multivalency', action="store_true")
    optional.add_argument('-nc','--no_constitutive', action="store_true", 
                    help='Exclude constitutive category from the output')
    optional.add_argument('-g',"--germsdir", type=str, default=os.getcwd(), nargs='?',
                        help='directory for where to find germs.R for multivalency analysis eg. /Users/Bellinda/repos/germs [DEFAULT current directory]')
    parser._action_groups.append(optional)
    args = parser.parse_args()

    logging.info(args)

    return(
        args.inputsplice,
        args.inputxlsites,
        args.genomefasta,
        args.fastaindex,
        args.outputpath,
        args.window,
        args.smoothing,
        args.minctrl,
        args.maxctrl,
        args.maxincl,
        args.maxfdr,
        args.maxenh,
        args.minsil,
        args.multivalency,
        args.germsdir,
        args.no_constitutive
        )

def df_apply(col_fn, *col_names):
    def inner_fn(df):
        cols = [df[col] for col in col_names]
        return col_fn(*cols)
    return inner_fn

def get_ss_bed(df, pos_col, neg_col):
    ss_pos = df.loc[df['strand'] == "+", ['chr', pos_col, pos_col, 'category', 'FDR', 'strand']]
    ss_pos.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_pos.start = ss_pos.start.transform(lambda x: x-1)

    ss_n = df.loc[df['strand'] == "-", ['chr', neg_col, neg_col, 'category', 'FDR', 'strand']]
    ss_n.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_n.end = ss_n.end.transform(lambda x: x+1)

    ss = pd.concat([ss_pos, ss_n])

    return ss

def get_coverage_plot(xl_bed, df, fai, window, exon_categories, label):
    """Return coverage of xl_bed items around df features extended by windows"""
    df = df.loc[df.name != "."]
    xl_bed = pbt.BedTool(xl_bed).sort()
    pbt_df = pbt.BedTool.from_dataframe(df[['chr', 'start', 'end', 'name', 'score', 'strand']]).sort().slop(l=window, r=window, s=True, g=fai)    
    df_coverage = pbt_df.coverage(b=xl_bed, **{'sorted': True, 's': True, 'd': True, 'nonamecheck': True}).to_dataframe()[['thickStart', 'thickEnd', 'strand', 'name']]
    df_coverage.rename(columns=({'thickStart':'position','thickEnd':'coverage'}), inplace=True)

    df_plot = df_coverage
    
    df_plot.loc[df_plot.strand=='+', 'position'] = df_plot['position'].astype('int32')
    df_plot.loc[df_plot.strand=='-', 'position'] = abs(2 * window + 2 - df_plot['position'])

    df_plot = df_plot.loc[df_plot.name != "."]
    df_plot = df_plot.groupby(['name','position'], as_index=False).agg({'coverage':'sum'})

    exon_cat = pd.DataFrame({'name':exon_categories.index, 'number_exons':exon_categories.values})
    df_plot = df_plot.merge(exon_cat, how="left")

    df_plot['norm_coverage'] = np.where(df_plot['coverage'] == 0, df_plot['coverage'], df_plot['coverage']/df_plot['number_exons'])

    # fisher test [covered exons, non-covered exon, control covered exons, control non-covered exons]
    df_plot_ctrl = df_plot.loc[df_plot.name == "control"][["position","coverage","number_exons"]]
    df_plot_ctrl.columns = ["position","control_coverage","control_number_exons"]
    df_plot = df_plot.merge(df_plot_ctrl, how="left")
    df_plot['control_norm_coverage'] = df_plot["control_coverage"] / df_plot["control_number_exons"]
    # give 0 values a small pseudo-count so that fold change can be calculated
    df_plot.loc[df_plot['control_norm_coverage'] == 0, ['control_norm_coverage']] = 0.000001
    df_plot['fold_change'] = df_plot["norm_coverage"] / df_plot["control_norm_coverage"]
    
    contingency_table = list(zip(df_plot['coverage'], df_plot['number_exons']-df_plot['coverage'], df_plot['control_coverage'], df_plot['control_number_exons'] - df_plot['control_coverage']))
    contingency_table = [ np.array(table).reshape(2,2) for table in contingency_table ]
    df_plot['pvalue'] = [ stats.fisher_exact(table)[1] for table in contingency_table ]

    df_plot['-log10pvalue'] = np.log10(1/df_plot['pvalue'])
    df_plot['label'] = label

    df_plot.loc[df_plot['fold_change'] < 1, ['-log10pvalue']] = df_plot['-log10pvalue'] * -1

    df_plot['-log10pvalue_smoothed'] = df_plot['-log10pvalue'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)

    return df_plot

def set_legend_text(legend, exon_categories, original_counts=None):
    """
    Set the legend text with optional subset information.
    
    Parameters:
    ----------
    legend : matplotlib.legend.Legend
        The legend object to modify
    exon_categories : pandas.Series
        Series containing category counts
    original_counts : dict, optional
        Dictionary with original counts before subsetting
    """
    # Convert exon_categories to a DataFrame for easier handling
    exon_cat = pd.DataFrame({'name': exon_categories.index, 'number_exons': exon_categories.values})
    
    # Determine which categories are present
    categories = exon_cat['name'].unique()
    legend.set_bbox_to_anchor([1.08, 0.75])
    legend.set_title("")
    
    # Track index in legend texts
    legend_idx = 0
    
    # Set text for each category in a specific order
    for category in ['constituitive', 'control', 'enhanced', 'silenced']:
        if category in categories:
            count = exon_cat[exon_cat['name'] == category]['number_exons'].values[0]
            text = f"{category.capitalize()} ({count}"
            
            # Add subset information if applicable
            if original_counts and category in original_counts:
                orig_count = original_counts[category]
                if orig_count > count:
                    text += f", subset from {orig_count}"
            
            text += ")"
            legend.texts[legend_idx].set_text(text)
            legend_idx += 1       

def add_enrichment_marker(fig, ax):
    """
    Add text-based enrichment/depletion marker with arrows on the right side of the figure,
    centered around y=0, dynamically accounting for the position of y=0.
    """
    # Get the y-axis limits of the current plot (ax)
    y_min, y_max = ax.get_ylim()
    
    # Position of y=0 relative to the y-axis
    y_zero_position = 0  # We're centering around y=0, no matter where it is on the y-axis
    
    # Calculate the height and position of the marker axis
    marker_height = 0.4  # Adjust this based on your plot, it defines how tall the marker is
    
    # Calculate the normalized position of the marker to be centered around y=0
    normalized_bottom = (y_zero_position - y_min) / (y_max - y_min) - marker_height / 2
    
    # Create a new axis on the right side of the figure
    marker_ax = fig.add_axes([0.94, normalized_bottom, 0.06, marker_height], label=f"enrichment_marker_{ax.get_title()}")
    
    # Remove all axis elements
    marker_ax.set_xticks([])
    marker_ax.set_yticks([])
    for spine in marker_ax.spines.values():
        spine.set_visible(False)
    
    # Position the arrows with text centered at the middle of the marker
    marker_ax.text(0.5, 0.5, "↑\nenriched\n\n\n↓\ndepleted",
                   ha='center', va='center', fontsize=11, 
                   transform=marker_ax.transAxes)
    
    return marker_ax



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

    mdf[['label','roname']] = mdf['sequence_name'].str.split('XX',expand=True)

    # GET MOST CONTRIBUTING KMERS
    # Add 'exon_type' column by splitting 'sequence_name' on 'XX'
    mdf['exon_type'] = mdf['sequence_name'].str.split('XX').str[0]
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

    mdf = mdf.groupby(['label','position'], as_index=False).agg({'smoothed_kmer_multivalency':'mean'})
    top_kmers_df = top_kmers_df.groupby(['label','position','exon_type','kmer'], as_index=False).agg({'smoothed_kmer_multivalency':'mean'})

    mdf['type'] = type
    top_kmers_df['type'] = type

    mdf = mdf.loc[(mdf.label != ".") & (pd.notnull(mdf.label)) & (mdf.label != "None")]
    top_kmers_df = top_kmers_df.loc[(top_kmers_df.label != ".") & (pd.notnull(top_kmers_df.label)) & (top_kmers_df.label != "None")]

    return mdf,top_kmers_df

def run_rna_map(de_file, xl_bed, genome_fasta, fai, window, smoothing, 
        min_ctrl, max_ctrl, max_inclusion, max_fdr, max_enh, min_sil, output_dir, multivalency, germsdir, no_constitutive
       #n_exons = 150, n_samples = 300, z_test=False
       ):
    name = de_file.split('/')[-1].replace('.txt', '').replace('.gz', '')
    df_fai = pd.read_csv(fai, sep='\t', header=None)
    chroms = set(df_fai[0].values)
    rmats = pd.read_csv(de_file, sep='\t')

    if 'exonStart_0base' in rmats.columns:
        rmats = rmats[rmats['chr'].isin(chroms)]
        rmats['inclusion'] = (rmats.IncLevel1.str.split(',') + rmats.IncLevel2.str.split(','))
        # code below for mean inclusion
        # rmats['inclusion'] = rmats['inclusion'].apply(lambda x: sum([float(y) for y in x if y != 'NA']) / len(x))
        # replaces with max inclusion
        rmats['inclusion'] = rmats['inclusion'].apply(lambda x: max([float(y) for y in x if y != 'NA']))
        df_rmats =  rmats.loc[ : ,['chr', 'exonStart_0base', 'exonEnd', 'FDR', 'IncLevelDifference', 'strand', 'inclusion', 
                                   'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']].rename(
            columns={'IncLevelDifference': 'dPSI', 'inclusion':'maxPSI'}).reset_index()
   
        
        # to deduplicate, first select the most extreme dPSI value for every exon (keep ties, they will be resolved by the hierarchy)
        df_rmats = df_rmats[df_rmats.groupby(['chr', 'exonStart_0base', 'exonEnd', 'strand'])['dPSI'].apply(lambda x: abs(x).rank(ascending=False) < 2)]


        # then apply hierarchy to decide which category exons belong to
        if no_constitutive:
            conditions = [
                (df_rmats["dPSI"].gt(min_sil) & df_rmats["FDR"].lt(max_fdr)), # silenced
                (df_rmats["dPSI"].lt(max_enh) & df_rmats["FDR"].lt(max_fdr)), # enhanced
                (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl))# control
            ]
            choices = ["silenced", "enhanced", "control"]
        else:
            conditions = [
                (df_rmats["dPSI"].gt(min_sil) & df_rmats["FDR"].lt(max_fdr)), # silenced
                (df_rmats["dPSI"].lt(max_enh) & df_rmats["FDR"].lt(max_fdr)), # enhanced
                (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl) & df_rmats["maxPSI"].gt(max_inclusion)), # constituitive
                (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl))# control
            ]
            choices = ["silenced", "enhanced", "constituitive", "control"]

        df_rmats["category"] = np.select(conditions, choices, default=None)

        df_rmats.to_csv(f'{output_dir}/{name}_RMATS_with_categories.tsv', sep="\t")

        exon_categories = df_rmats.groupby('category').size()
        #exon_categories.columns = ["name", "exon_number"]
        logging.info("Exons in each category:")
        logging.info(exon_categories)
        logging.info("Total categorised deduplicated exons: " + str(df_rmats.shape[0]))

        ### Some warning messages ###
        if exon_categories.loc["control"] == 0:
            logging.info("Warning! There are no control exons. Try changing thresholds or input file and run again.")
            sys.exit()
        
        if exon_categories.loc["enhanced"] == 0 and exon_categories.loc["silenced"] == 0:
            logging.info('Warning! There are no regulated exons, try changing filtering parameters or file and run again.')
            sys.exit()
        
        # Apply subsetting for the plotting
        # Get the count of each category
        category_counts = df_rmats['category'].value_counts()
        # Store original counts
        original_counts = {cat: count for cat, count in category_counts.items()}
    
        # Find the maximum count of enhanced or silenced
        target_count = 0
        if 'enhanced' in category_counts and 'silenced' in category_counts:
            target_count = max(category_counts['enhanced'], category_counts['silenced'])
        elif 'enhanced' in category_counts:
            target_count = category_counts['enhanced']
        elif 'silenced' in category_counts:
            target_count = category_counts['silenced']
    
        # Subset control exons
        if 'control' in category_counts and category_counts['control'] > target_count > 0:
            control_indices = df_rmats[df_rmats['category'] == 'control'].index
            # Randomly select indices to keep
            control_indices_to_keep = np.random.choice(control_indices, target_count, replace=False)
            # Create a mask for rows to drop
            drop_mask = df_rmats.index.isin(control_indices) & ~df_rmats.index.isin(control_indices_to_keep)
            # Drop the rows
            df_rmats = df_rmats[~drop_mask]
            logging.info(f"Randomly subsetted control exons from {category_counts['control']} to {target_count}")
    
        # Subset constitutive exons if they exist and not excluded
        if not no_constitutive and 'constituitive' in category_counts and category_counts['constituitive'] > target_count > 0:
            const_indices = df_rmats[df_rmats['category'] == 'constituitive'].index
            # Randomly select indices to keep
            const_indices_to_keep = np.random.choice(const_indices, target_count, replace=False)
            # Create a mask for rows to drop
            drop_mask = df_rmats.index.isin(const_indices) & ~df_rmats.index.isin(const_indices_to_keep)
            # Drop the rows
            df_rmats = df_rmats[~drop_mask]
            logging.info(f"Randomly subsetted constitutive exons from {category_counts['constituitive']} to {target_count}")

        exon_categories = df_rmats.groupby('category').size()

        ####### Exon lengths #######
        df_rmats["regulated_exon_length"] = df_rmats['exonEnd'] - df_rmats['exonStart_0base']
        df_rmats["upstream_exon_length"] = df_rmats['upstreamEE'] - df_rmats['upstreamES']
        df_rmats["downstream_exon_length"] = df_rmats['downstreamEE'] - df_rmats['downstreamES']

        exon_length_df = df_rmats[["regulated_exon_length","upstream_exon_length","downstream_exon_length","category"]]

        exon_length_df = exon_length_df .melt(id_vars=["category"], var_name="exon_type", value_name="exon_length")
        
        palette_exon_len = [colors_dict['ctrl'], colors_dict['const'], colors_dict['enh'], 
                            colors_dict['enhrest'], colors_dict['sil'], colors_dict['silrest'], colors_dict['all']]

        sns.set(rc={'figure.figsize':(15, 5)})
        sns.set_style("whitegrid")
        g = sns.catplot(data=exon_length_df, x='category', y='exon_length',col='exon_type', 
                    kind='box', col_wrap=3, showfliers=False,
                    col_order=["upstream_exon_length","regulated_exon_length","downstream_exon_length"],
                    order=["control","constituitive","enhanced","enhanced_rest","silenced","silenced_rest"],
                    palette = palette_exon_len)
        titles = ["Upstream Exon", "Middle Exon", "Downstream exon"]
        for ax, title in zip(g.axes.flat, titles):
            ax.set_title(title)
        g.set_xticklabels(rotation=45)
        g.set(xlabel=None)
        g.axes[0].set_ylabel('Exon length (bp)')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{name}_exon_length.pdf')
        pbt.helpers.cleanup()


        ### The coverage plot ###

        middle_3ss_bed = get_ss_bed(df_rmats,'exonStart_0base','exonEnd')
        middle_5ss_bed = get_ss_bed(df_rmats,'exonEnd','exonStart_0base')
        downstream_3ss_bed = get_ss_bed(df_rmats,'downstreamES','downstreamEE')
        downstream_5ss_bed = get_ss_bed(df_rmats,'downstreamEE','downstreamES')
        upstream_3ss_bed = get_ss_bed(df_rmats,'upstreamES','upstreamEE')
        upstream_5ss_bed = get_ss_bed(df_rmats,'upstreamEE','upstreamES')

        if xl_bed is not None:
            middle_3ss = get_coverage_plot(xl_bed, middle_3ss_bed, fai, window, exon_categories, 'middle_3ss')
            middle_5ss = get_coverage_plot(xl_bed, middle_5ss_bed, fai, window, exon_categories, 'middle_5ss')
            downstream_3ss = get_coverage_plot(xl_bed, downstream_3ss_bed, fai, window, exon_categories, 'downstream_3ss')
            downstream_5ss = get_coverage_plot(xl_bed, downstream_5ss_bed, fai, window, exon_categories, 'downstream_5ss')
            upstream_3ss = get_coverage_plot(xl_bed, upstream_3ss_bed, fai, window, exon_categories, 'upstream_3ss')
            upstream_5ss = get_coverage_plot(xl_bed, upstream_5ss_bed, fai, window, exon_categories, 'upstream_5ss')


            plotting_df = pd.concat([middle_3ss, middle_5ss, downstream_3ss, downstream_5ss, upstream_3ss, upstream_5ss])
            plotting_df.to_csv(f'{output_dir}/{name}_RNAmap.tsv', sep="\t")

            #sns.set(rc={'figure.figsize':(7, 5)})
            sns.set_style("whitegrid")

            g = sns.relplot(data=plotting_df, x='position', y='-log10pvalue_smoothed', hue='name', col='label', facet_kws={"sharex":False},
                        kind='line', col_wrap=6, height=5, aspect=4/5,
                        col_order=["upstream_3ss","upstream_5ss","middle_3ss","middle_5ss","downstream_3ss","downstream_5ss"])
            titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS", "Downstream 5'SS"]
            for ax, title in zip(g.axes.flat, titles):
                ax.set_title(title)
                ax.axhline(y=0, color='k', alpha=0.2, linewidth=0.5)
                fig = plt.gcf()
                marker_ax = add_enrichment_marker(fig, ax) 

            g.set(xlabel='')
            g.axes[0].set_ylabel('-log10(p value) enrichment / control')

            sns.move_legend(
                g, "upper right",  # Position: 'upper left' relative to bbox
                bbox_to_anchor=(1, 2),  # Outside the plot area to the right
                ncol=1,  # Number of columns in the legend
                title=None,  # Set title if needed
                frameon=False  # Remove frame around the legend
            )
            leg = g._legend
            set_legend_text(leg, exon_categories, original_counts)

			
            ### Add exon-intron drawing below line plots ###
            # for calculating rectangle size
            rect_fraction = 1 / ((window + 50) / 50)

            ax = g.axes[0]
            ax.set_xlim([0, window+50])
            a=ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(a))
            a = np.arange(0-window, 51, 50)
            a = list(map(str, a))
            a[0] = ""
            a[-1] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            ax.set_xlim([window-50, window*2])
            a=ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(a))
            a = np.arange(-50,window + 1, 50)
            a = list(map(str, a))
            a[0] = ""
            a[-1] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)



            ax = g.axes[1]
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[2]
            ax.set_xlim([0, window+50])
            a=ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(a))
            a = np.arange(0-window, 51, 50)
            a = list(map(str, a))
            a[0] = ""
            a[-1] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[2]
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[3]
            ax.set_xlim([window-50, window*2])
            a=ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(a))
            a = np.arange(-50,window + 1, 50)
            a = list(map(str, a))
            a[0] = ""
            a[-1] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[3]
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[4]
            ax.set_xlim([0, window+50])
            a=ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(a))
            a = np.arange(0-window, 51, 50)
            a = list(map(str, a))
            a[0] = ""
            a[-1] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[5]
            ax.set_xlim([window-50, window*2])
            a=ax.get_xticks().tolist()
            ax.xaxis.set_major_locator(mticker.FixedLocator(a))
            a = np.arange(-50,window + 1, 50)
            a = list(map(str, a))
            a[0] = ""
            a[-1] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[5]
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)


            # Adjust subplot spacing and save with enough room for legend and arrows
            plt.subplots_adjust(wspace=0.05)  
            plt.savefig(f'{output_dir}/{name}_RNAmap_-log10pvalue.pdf', 
                    bbox_extra_artists=([leg, rect, marker_ax]),
                    bbox_inches='tight',
                    pad_inches=0.8)
            pbt.helpers.cleanup()

        if multivalency:
            ### Get multivalency scores ###
            rect_fraction = 1 / ((window + 50) / 50)

            middle_3ss_mdf = get_multivalency_scores(middle_3ss_bed, fai, window, genome_fasta, output_dir, name, 'middle_3ss',germsdir)[0]
            middle_5ss_mdf = get_multivalency_scores(middle_5ss_bed, fai, window, genome_fasta, output_dir, name, 'middle_5ss',germsdir)[0]
            downstream_3ss_mdf = get_multivalency_scores(downstream_3ss_bed, fai, window, genome_fasta, output_dir, name, 'downstream_3ss',germsdir)[0]
            downstream_5ss_mdf = get_multivalency_scores(downstream_5ss_bed, fai, window, genome_fasta, output_dir, name, 'downstream_5ss',germsdir)[0]
            upstream_3ss_mdf = get_multivalency_scores(upstream_3ss_bed, fai, window, genome_fasta, output_dir, name, 'upstream_3ss',germsdir)[0]
            upstream_5ss_mdf = get_multivalency_scores(upstream_5ss_bed, fai, window, genome_fasta, output_dir, name, 'upstream_5ss',germsdir)[0]

            plotting_df = pd.concat([middle_3ss_mdf, middle_5ss_mdf, downstream_3ss_mdf, downstream_5ss_mdf, upstream_3ss_mdf, upstream_5ss_mdf])
            plotting_df.to_csv(f'{output_dir}/{name}_RNAmap_multivalency.tsv', sep="\t")

            plt.figure()
            sns.set_style("whitegrid")
            g = sns.relplot(data=plotting_df, x='position', y='smoothed_kmer_multivalency', hue='label', col='type', facet_kws={"sharex":False},
                        kind='line', col_wrap=6, height=5, aspect=3.5/5,
                        col_order=["upstream_3ss","upstream_5ss","middle_3ss","middle_5ss","downstream_3ss","downstream_5ss"])
            titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS", "Downstream 5'SS"]
            for ax, title in zip(g.axes.flat, titles):
                ax.set_title(title)
            g.set(xlabel='')
            g.axes[0].set_ylabel('mean smoothed kmer multivalency')
            leg = g._legend
            set_legend_text(leg, exon_categories, original_counts)

            ax = g.axes[0]
            ax.set_xlim([window, (2*window)+50])
            ax.set_ylim(ymin=1)
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[2]
            ax.set_xlim([window, (2*window)+50])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[3]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[4]
            ax.set_xlim([window, (2*window)+50])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[5]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            plt.subplots_adjust(wspace=0.01)
            plt.savefig(f'{output_dir}/{name}_RNAmap_multivalency.pdf',
                bbox_extra_artists=([leg,rect]),
                bbox_inches='tight',
                pad_inches=0.5)
            pbt.helpers.cleanup()



        # DO THE SAME NOW FOR THE KMERS
            middle_3ss_mdf = get_multivalency_scores(middle_3ss_bed, fai, window, genome_fasta, output_dir, name, 'middle_3ss',germsdir)[1]
            middle_5ss_mdf = get_multivalency_scores(middle_5ss_bed, fai, window, genome_fasta, output_dir, name, 'middle_5ss',germsdir)[1]
            downstream_3ss_mdf = get_multivalency_scores(downstream_3ss_bed, fai, window, genome_fasta, output_dir, name, 'downstream_3ss',germsdir)[1]
            downstream_5ss_mdf = get_multivalency_scores(downstream_5ss_bed, fai, window, genome_fasta, output_dir, name, 'downstream_5ss',germsdir)[1]
            upstream_3ss_mdf = get_multivalency_scores(upstream_3ss_bed, fai, window, genome_fasta, output_dir, name, 'upstream_3ss',germsdir)[1]
            upstream_5ss_mdf = get_multivalency_scores(upstream_5ss_bed, fai, window, genome_fasta, output_dir, name, 'upstream_5ss',germsdir)[1]

            plotting_df = pd.concat([middle_3ss_mdf, middle_5ss_mdf, downstream_3ss_mdf, downstream_5ss_mdf, upstream_3ss_mdf, upstream_5ss_mdf])
            plotting_df.to_csv(f'{output_dir}/{name}_RNAmap_TOP10KMER_multivalency.tsv', sep="\t")

            plt.figure()
            sns.set_style("whitegrid")
            g = sns.relplot(data=plotting_df[plotting_df['exon_type'] == 'silenced'], x='position', y='smoothed_kmer_multivalency', hue='kmer', col='type', facet_kws={"sharex":False},
                        kind='line', col_wrap=6, height=5, aspect=3.5/5,
                        col_order=["upstream_3ss","upstream_5ss","middle_3ss","middle_5ss","downstream_3ss","downstream_5ss"])
            titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS", "Downstream 5'SS"]
            for ax, title in zip(g.axes.flat, titles):
                ax.set_title(title)
            g.set(xlabel='')
            g.axes[0].set_ylabel('mean smoothed kmer multivalency')
            leg = g._legend
            set_legend_text(leg, exon_categories, original_counts)

            ax = g.axes[0]
            ax.set_xlim([window, (2*window)+50])
            ax.set_ylim(ymin=1)
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[2]
            ax.set_xlim([window, (2*window)+50])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[3]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[4]
            ax.set_xlim([window, (2*window)+50])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[5]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            plt.subplots_adjust(wspace=0.01)
            plt.savefig(f'{output_dir}/{name}_RNAmap_silencedKMER_multivalency.pdf',
                bbox_extra_artists=([leg,rect]),
                bbox_inches='tight',
                pad_inches=0.5)
            pbt.helpers.cleanup()

            plt.figure()
            sns.set_style("whitegrid")
            g = sns.relplot(data=plotting_df[plotting_df['exon_type'] == 'enhanced'], x='position', y='smoothed_kmer_multivalency', hue='kmer', col='type', facet_kws={"sharex":False},
                        kind='line', col_wrap=6, height=5, aspect=3.5/5,
                        col_order=["upstream_3ss","upstream_5ss","middle_3ss","middle_5ss","downstream_3ss","downstream_5ss"])
            titles = ["Upstream 3'SS", "Upstream 5'SS", "Middle 3'SS", "Middle 5'SS", "Downstream 3'SS", "Downstream 5'SS"]
            for ax, title in zip(g.axes.flat, titles):
                ax.set_title(title)
            g.set(xlabel='')
            g.axes[0].set_ylabel('mean smoothed kmer multivalency')
            leg = g._legend
            set_legend_text(leg, exon_categories, original_counts)
            ax = g.axes[0]
            ax.set_xlim([window, (2*window)+50])
            ax.set_ylim(ymin=1)
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[1]
            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[2]
            ax.set_xlim([window, (2*window)+50])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[3]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="midnightblue", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[4]
            ax.set_xlim([window, (2*window)+50])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-window,50, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)

            rect = matplotlib.patches.Rectangle(
                xy=(1 - rect_fraction, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            ax = g.axes[5]
            ax.set_xlim([(2*window)-50, 3*window])
            start, end = ax.get_xlim()
            ax.xaxis.set_ticks(np.arange(start, end, 50))
            a=ax.get_xticks().tolist()
            a = np.arange(-50,window, 50)
            a = list(map(str, a))
            a[0] = ""
            ax.set_xticklabels(a)
            rect = matplotlib.patches.Rectangle(
                xy=(0, -0.2), width=rect_fraction, height=.1,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            rect = matplotlib.patches.Rectangle(
                xy=(rect_fraction, -0.15), width=1 - rect_fraction, height=.001,
                color="slategrey", alpha=1,
                transform=ax.transAxes, clip_on=False,
                )
            ax.add_artist(rect)

            plt.subplots_adjust(wspace=0.01)
            plt.savefig(f'{output_dir}/{name}_RNAmap_enhancedKMER_multivalency.pdf',
                bbox_extra_artists=([leg,rect]),
                bbox_inches='tight',
                pad_inches=0.5)
            pbt.helpers.cleanup()
        sys.exit()


 
if __name__=='__main__':
    (
        de_file,
        xl_bed,
        genome_fasta,
        fai,
        output_folder,
        window,
        smoothing, 
        min_ctrl,
        max_ctrl,
        max_inclusion,
        max_fdr,
        max_enh,
        min_sil,
        multivalency,
        germsdir,
        no_constitutive
    ) = cli()
    
    log_filename, start_time, logger = setup_logging(output_folder)
    logging.info(f"Log file created: {log_filename}")

    try:
        run_rna_map(de_file, xl_bed, genome_fasta, fai, window, smoothing, 
            min_ctrl, max_ctrl, max_inclusion, max_fdr, max_enh, min_sil, output_folder, multivalency, germsdir, no_constitutive)
    finally:
        # Log runtime at the end
        log_runtime(start_time, logger)
        
        # Explicitly flush logs before shutting down
        for handler in logger.handlers:
            handler.flush()

        # Now safe to shut down logging
        logging.shutdown()
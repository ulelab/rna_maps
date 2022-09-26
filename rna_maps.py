import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import pandas as pd
import numpy as np
import pybedtools as pbt
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys
import os
import argparse

sns.set_style("whitegrid")
colors_dict = {'all': '#D3D3D3', 'ctrl': '#408F76', 'enh': '#F30C08', 'sil': '#005299', 'enhrest': '#FFB122', 'silrest': '#6DC2F5', 'const': '#666666'}
linewidth = 3
dashes = False

def cli():
    parser = argparse.ArgumentParser(description='Plot CLIP crosslinks around regulated exons to study position-dependent impact on pre-mRNA splicing.')
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i',"--inputsplice", type=str, required=True,
                        help='quantification of differential splicing produced by rMATS')
    required.add_argument('-x',"--inputxlsites", type=str, required=True,
                        help='CLIP crosslinks in BED file format')
    required.add_argument('-f',"--fastaindex", type=str, required=True,
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
    parser._action_groups.append(optional)
    args = parser.parse_args()
    print(args)

    return(
        args.inputsplice,
        args.inputxlsites,
        args.fastaindex,
        args.outputpath,
        args.window,
        args.smoothing,
        args.minctrl,
        args.maxctrl,
        args.maxincl,
        args.maxfdr,
        args.maxenh,
        args.minsil
        )

def get_coverage_plot(xl_bed, df, fai, window):
    """Return coverage of xl_bed items around df features extended by windows"""
    xl_bed = pbt.BedTool(xl_bed).sort()
    pbt_df = pbt.BedTool.from_dataframe(df[['chr', 'start', 'end', 'name', 'score', 'strand']]).sort().slop(l=window, r=window, s=True, g=fai)    
    df_coverage = pbt_df.coverage(b=xl_bed, **{'sorted': True, 's': True, 'd': True, 'nonamecheck': True}).to_dataframe()[
        ['thickStart', 'thickEnd', 'strand']].rename(columns={'thickStart': 'position', 'thickEnd': 'coverage'})
    df_plot = df_coverage.groupby(['position', 'strand']).sum().reset_index()
    df_plot.loc[df_plot.strand=='+', 'map'] = df_plot['position'].astype('int32')
    df_plot.loc[df_plot.strand=='-', 'map'] = abs(2 * window + 2 - df_plot['position'])
    df_plot.map = df_plot.map.astype('int32')
    df_plot = df_plot[['coverage', 'map']].groupby('map').sum()
    df_plot_norm = df_plot / len(df)
    return df_plot_norm, df_coverage, df_plot
            
            
def get_exon_dist(len_df_in, df_coverage, window):
    p = 2 * window + 1
    no_xl_per_exon = {}
    for n in range(len_df_in):
        no_xl = df_coverage.loc[n * p: n * p + p, 'coverage'].sum()
        no_xl_per_exon[n] = no_xl
    return no_xl_per_exon


def get_3ss5ss_exons(df_exons):
    exons_stranded = df_exons.copy()
    exons_stranded.loc[df_exons[5] == '-', 1] = df_exons[2] # replace start with end
    exons_stranded.loc[df_exons[5] == '-', 2] = df_exons[1] # replace end with start
    exons_stranded.loc[df_exons[5] == '-', 7] = df_exons[8] # replace start with end upstream exon
    exons_stranded.loc[df_exons[5] == '-', 8] = df_exons[7] # replace end with start upstream exon
    exons_stranded.loc[df_exons[5] == '-', 9] = df_exons[10] # replace start with end downstream exon
    exons_stranded.loc[df_exons[5] == '-', 10] = df_exons[9] # replace end with start downstream exon
    exons_3ss = exons_stranded[[0, 1, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14]].copy() # takes all the starts for 3ss
    exons_5ss = exons_stranded[[0, 2, 3, 4, 5, 6, 8, 10, 11, 12, 13, 14]].copy() # takes all the ends for 5ss
    exons_3ss[2] = exons_3ss[1] + 1 # calculates end from start
    exons_5ss[1] = exons_5ss[2] # assigns start to end for 5ss
    exons_5ss[2] = exons_5ss[1] + 1 # calulates end from start
    exons_3ss[8] = exons_3ss[7] + 1 # calculates end from start for upstream exons
    exons_5ss[7] = exons_5ss[8] # assigns start to end for 5ss for upstream exons
    exons_5ss[8] = exons_5ss[7] + 1 # calulates end from start for upstream exons
    exons_3ss[10] = exons_3ss[9] + 1 # calculates end from start for downstream exons
    exons_5ss[9] = exons_5ss[10] # assigns start to end for 5ss for downstream exons
    exons_5ss[10] = exons_5ss[9] + 1 # calulates end from start for downstream exons
    return exons_3ss[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]], exons_5ss[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]]




def run_rna_map(de_file, xl_bed, fai, window, smoothing, 
        min_ctrl, max_ctrl, max_inclusion, max_fdr, max_enh, min_sil, output_dir,
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
        conditions = [
            (df_rmats["dPSI"].gt(min_sil) & df_rmats["FDR"].lt(max_fdr)), # silenced
            (df_rmats["dPSI"].lt(max_enh) & df_rmats["FDR"].lt(max_fdr)), # enhanced
            (df_rmats["dPSI"].gt(min_sil)), # silenced rest
            (df_rmats["dPSI"].lt(max_enh)), # enhanced rest
            (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl) & df_rmats["maxPSI"].gt(max_inclusion)), # constituitive
            (df_rmats["dPSI"].gt(min_ctrl) & df_rmats["dPSI"].lt(max_ctrl))# control
        ]

        choices = ["silenced", "enhanced", "silenced_rest", "enhanced_rest", "constituitive", "control"]

        df_rmats["category"] = np.select(conditions, choices, default=None)

        # logging info
        print("Exons in each category:")
        print(df_rmats.groupby('category').size())
        print("Total categorised deduplicated exons: ", str(df_rmats.shape[0]))

        ####### Exon lengths #######
        df_rmats["exon_length"]
        df_rmats["upstream_exon_length"]
        df_rmats["downstream_exon_length"]

        sys.exit

    if de_source == 'rmats':
        df_rmats_enh_3ss, df_rmats_enh_5ss = get_3ss5ss_exons(df_rmats_enh)
        print(f'There are: {len(df_rmats_enh_3ss)} 3ss and {len(df_rmats_enh_5ss)} 5ss enhanced exons')
        df_rmats_sil_3ss, df_rmats_sil_5ss = get_3ss5ss_exons(df_rmats_sil)
        print(f'There are: {len(df_rmats_sil_3ss)} 3ss and {len(df_rmats_sil_5ss)} 5ss silenced exons')
        df_rmats_enhrest_3ss, df_rmats_enhrest_5ss = get_3ss5ss_exons(df_rmats_enhrest)
        print(f'There are: {len(df_rmats_enhrest_3ss)} 3ss and {len(df_rmats_enhrest_5ss)} 5ss possibly enhanced exons with FDR that is to hight')
        df_rmats_silrest_3ss, df_rmats_silrest_5ss = get_3ss5ss_exons(df_rmats_silrest)
        print(f'There are: {len(df_rmats_silrest_3ss)} 3ss and {len(df_rmats_silrest_5ss)} 5ss possibly silenced exons with FDR that is to hight')
        df_rmats_ctrl_3ss, df_rmats_ctrl_5ss = get_3ss5ss_exons(df_rmats_ctrl)
        print(f'There are: {len(df_rmats_ctrl_3ss)} 3ss and {len(df_rmats_ctrl_5ss)} 5ss control exons')
        df_rmats_const_3ss, df_rmats_const_5ss = get_3ss5ss_exons(df_rmats_const)
        print(f'There are: {len(df_rmats_const_3ss)} 3ss and {len(df_rmats_const_5ss)} 5ss constitutive exons')
        
        col_rename = {0: 'chr', 1: 'start', 2: 'end', 3: 'name', 4: 'score', 5: 'strand', 6:'inclusion', 7: 'upstream_es', 
            8: 'upstream_ee', 9: 'downstream_es', 10: 'downstream_ee' , 11: 'exon_len', 12: 'upstream_exon_len', 13: 'downstream_exon_len'
                     }
        
        df_rmats_enh_3ss = df_rmats_enh_3ss.rename(columns=col_rename)
        df_rmats_enh_5ss = df_rmats_enh_5ss.rename(columns=col_rename)
        df_rmats_sil_3ss = df_rmats_sil_3ss.rename(columns=col_rename)
        df_rmats_sil_5ss = df_rmats_sil_5ss.rename(columns=col_rename)
        df_rmats_enhrest_3ss = df_rmats_enhrest_3ss.rename(columns=col_rename)
        df_rmats_enhrest_5ss = df_rmats_enhrest_5ss.rename(columns=col_rename)
        df_rmats_silrest_3ss = df_rmats_silrest_3ss.rename(columns=col_rename)
        df_rmats_silrest_5ss = df_rmats_silrest_5ss.rename(columns=col_rename)
        df_rmats_ctrl_3ss = df_rmats_ctrl_3ss.rename(columns=col_rename)
        df_rmats_ctrl_5ss = df_rmats_ctrl_5ss.rename(columns=col_rename)
        df_rmats_const_3ss = df_rmats_const_3ss.rename(columns=col_rename)
        df_rmats_const_5ss = df_rmats_const_5ss.rename(columns=col_rename)

        index_selected = df_rmats_ctrl_3ss.index.union(df_rmats_enh_3ss.index.union(
            df_rmats_sil_3ss.index.union(df_rmats_enhrest_3ss.index.union(df_rmats_silrest_3ss.index.union(df_rmats_const_3ss.index)))))

                
    if len(df_rmats_ctrl_3ss) == 0:
        print('No control exons, try changing filtering parameters or file')
        return
    if len(df_rmats_enh_3ss) == 0 and len(df_rmats_sil_3ss) == 0:
        print('No regulated exons, try changing filtering parameters or file')
        return
#     if de_source == 'rmats':
#         exon_len_all = df_rmats[11]
#     elif de_source == 'whippet':
#         exon_len_all = df_rmats['exon_len']
    
    df_rmats_ctrl_3ss['exon_len'] = df_rmats_ctrl_3ss['end'] - df_rmats_ctrl_3ss['start']
    df_rmats_enh_3ss['exon_len'] = df_rmats_enh_3ss['end'] - df_rmats_enh_3ss['start']
    df_rmats_sil_3ss['exon_len'] = df_rmats_sil_3ss['end'] - df_rmats_sil_3ss['start']
    df_rmats_enhrest_3ss['exon_len'] = df_rmats_enhrest_3ss['end'] - df_rmats_enhrest_3ss['start']
    df_rmats_silrest_3ss['exon_len'] = df_rmats_silrest_3ss['end'] - df_rmats_silrest_3ss['start']
    df_rmats_const_3ss['exon_len'] = df_rmats_const_3ss['end'] - df_rmats_const_3ss['start']
    
    
    exon_len_all = list(df_rmats_ctrl_3ss['exon_len']) + list(df_rmats_enh_3ss['exon_len'] )+ list(df_rmats_sil_3ss['exon_len'])  \
        + list(df_rmats_enhrest_3ss['exon_len']) + list(df_rmats_silrest_3ss['exon_len']) + list(df_rmats_const_3ss['exon_len'])
    exon_len_ctrl = df_rmats_ctrl_3ss['exon_len']
    exon_len_enh = df_rmats_enh_3ss['exon_len']
    exon_len_sil = df_rmats_sil_3ss['exon_len']
    exon_len_enhrest = df_rmats_enhrest_3ss['exon_len']
    exon_len_silrest = df_rmats_silrest_3ss['exon_len']
    exon_len_const = df_rmats_const_3ss['exon_len']
    df_exon_len = pd.DataFrame(data={'enhanced': exon_len_enh, 'enhanced_rest': exon_len_enhrest, 'silenced': exon_len_sil, 
            'silenced_rest': exon_len_silrest, 'control': exon_len_ctrl, 
            'constitutive': exon_len_const, 'all_exons': exon_len_all})
    
    if de_source == 'rmats':
        upstream_exon_len_all = list(df_rmats_ctrl_3ss['upstream_exon_len']) + list(df_rmats_enh_3ss['upstream_exon_len']) + \
            list(df_rmats_sil_3ss['upstream_exon_len']) + list(df_rmats_enhrest_3ss['upstream_exon_len']) + \
            list(df_rmats_silrest_3ss['upstream_exon_len']) + list(df_rmats_const_3ss['upstream_exon_len'])
        upstream_exon_len_ctrl = df_rmats_ctrl_3ss['upstream_exon_len']
        upstream_exon_len_enh = df_rmats_enh_3ss['upstream_exon_len']
        upstream_exon_len_sil = df_rmats_sil_3ss['upstream_exon_len']
        upstream_exon_len_enhrest = df_rmats_enhrest_3ss['upstream_exon_len']
        upstream_exon_len_silrest = df_rmats_silrest_3ss['upstream_exon_len']
        upstream_exon_len_const = df_rmats_const_3ss['upstream_exon_len']
        df_upstream_exon_len = pd.DataFrame(data={'enhanced': upstream_exon_len_enh, 'enhanced_rest': upstream_exon_len_enhrest, 'silenced': upstream_exon_len_sil, 
            'silenced_rest': upstream_exon_len_silrest, 'control': upstream_exon_len_ctrl, 
            'constitutive': upstream_exon_len_const, 'all_exons': upstream_exon_len_all})

        downstream_exon_len_all = list(df_rmats_ctrl_3ss['downstream_exon_len']) + list(df_rmats_enh_3ss['downstream_exon_len']) + \
            list(df_rmats_sil_3ss['downstream_exon_len']) + list(df_rmats_enhrest_3ss['downstream_exon_len']) + \
            list(df_rmats_silrest_3ss['downstream_exon_len']) + list(df_rmats_const_3ss['downstream_exon_len'])
        downstream_exon_len_ctrl = df_rmats_ctrl_3ss['downstream_exon_len']
        downstream_exon_len_enh = df_rmats_enh_3ss['downstream_exon_len']
        downstream_exon_len_sil = df_rmats_sil_3ss['downstream_exon_len']
        downstream_exon_len_enhrest = df_rmats_enhrest_3ss['downstream_exon_len']
        downstream_exon_len_silrest = df_rmats_silrest_3ss['downstream_exon_len']
        downstream_exon_len_const = df_rmats_const_3ss['downstream_exon_len']
        df_downstream_exon_len = pd.DataFrame(data={  
             'enhanced': downstream_exon_len_enh, 'enhanced_rest': downstream_exon_len_enhrest, 'silenced': downstream_exon_len_sil, 
            'silenced_rest': downstream_exon_len_silrest, 'control': downstream_exon_len_ctrl, 
            'constitutive': downstream_exon_len_const, 'all_exons': downstream_exon_len_all})
        
    
    if de_source == 'rmats':
        sns.set(rc={'figure.figsize':(15, 5)})
        sns.set_style("whitegrid")
        palette_exon_len = [colors_dict['enh'], colors_dict['enhrest'], colors_dict['sil'], colors_dict['silrest'],
                            colors_dict['ctrl'], colors_dict['const'], colors_dict['all']]
        fig05, axs = plt.subplots(1, 3, sharey='row')
        ax0 = sns.boxplot(data=df_upstream_exon_len, palette=palette_exon_len,
            showfliers=False, ax=axs[0])
        ax0.set_xticklabels(ax0.get_xticklabels(), rotation=45)
        ax1 = sns.boxplot(data=df_exon_len, palette=palette_exon_len,
            showfliers=False, ax=axs[1])
        ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
        ax2 = sns.boxplot(data=df_downstream_exon_len, palette=palette_exon_len,
            showfliers=False, ax=axs[2])
        ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
        axs[0].set_title("Upstream exon length")
        axs[2].set_title("Downstream exon length")
        axs[0].set_ylabel('Exon length (bases)')
        plt.tight_layout()
        fig05.savefig(f'{output_dir}/{name}_exon_len.pdf')
    
    # if output_exons:
    #     name = de_file.split('/')[-1].split('_')[0]
    #     df_rmats_ctrl_3ss.to_csv(f'{output_dir}/{name}_{de_source}_ctrl_3ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_ctrl_5ss.to_csv(f'{output_dir}/{name}_{de_source}_ctrl_5ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_enh_3ss.to_csv(f'{output_dir}/{name}_{de_source}_enh_3ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_enh_5ss.to_csv(f'{output_dir}/{name}_{de_source}_enh_5ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_sil_3ss.to_csv(f'{output_dir}/{name}_{de_source}_sil_3ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_sil_5ss.to_csv(f'{output_dir}/{name}_{de_source}_sil_5ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_enhrest_3ss.to_csv(f'{output_dir}/{name}_{de_source}_enhrest_3ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_enhrest_5ss.to_csv(f'{output_dir}/{name}_{de_source}_enhrest_5ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_silrest_3ss.to_csv(f'{output_dir}/{name}_{de_source}_silrest_3ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_silrest_5ss.to_csv(f'{output_dir}/{name}_{de_source}_silrest_5ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_const_3ss.to_csv(f'{output_dir}/{name}_{de_source}_const_3ss.tsv', sep='\t', index=None, header=None)
    #     df_rmats_const_5ss.to_csv(f'{output_dir}/{name}_{de_source}_const_5ss.tsv', sep='\t', index=None, header=None)
    #     print('Exons saved in tsv files')
    

    col_bed = ['chr', 'start', 'end', 'name', 'score', 'strand']

    
    df_enh_3ss, df_coverage_enh_3ss, df_raw_enh_3ss = get_coverage_plot(xl_bed, df_rmats_enh_3ss[col_bed], fai, window)
    df_enh_3ss = df_enh_3ss.rename(columns={'coverage': 'enhanced'})
    df_enh_5ss, df_coverage_enh_5ss, df_raw_enh_5ss = get_coverage_plot(xl_bed, df_rmats_enh_5ss[col_bed], fai, window)
    df_enh_5ss = df_enh_5ss.rename(columns={'coverage': 'enhanced'})

    df_sil_3ss, df_coverage_sil_3ss, df_raw_sil_3ss = get_coverage_plot(xl_bed, df_rmats_sil_3ss[col_bed], fai, window)
    df_sil_3ss = df_sil_3ss.rename(columns={'coverage': 'silenced'})
    df_sil_5ss, df_coverage_sil_5ss, df_raw_sil_5ss = get_coverage_plot(xl_bed, df_rmats_sil_5ss[col_bed], fai, window)
    df_sil_5ss = df_sil_5ss.rename(columns={'coverage': 'silenced'})
    
    df_ctrl_3ss, df_coverage_ctrl_3ss, df_raw_ctrl_3ss = get_coverage_plot(xl_bed, df_rmats_ctrl_3ss[col_bed], fai, window)
    df_ctrl_3ss = df_ctrl_3ss.rename(columns={'coverage': 'control'})
    df_ctrl_5ss, df_coverage_ctrl_5ss, df_raw_ctrl_5ss = get_coverage_plot(xl_bed, df_rmats_ctrl_5ss[col_bed], fai, window)
    df_ctrl_5ss = df_ctrl_5ss.rename(columns={'coverage': 'control'})
    
    df_const_3ss, df_coverage_const_3ss, df_raw_const_3ss = get_coverage_plot(xl_bed, df_rmats_const_3ss[col_bed], fai, window)
    df_const_3ss = df_const_3ss.rename(columns={'coverage': 'const'})
    df_const_5ss, df_coverage_const_5ss, df_raw_const_5ss = get_coverage_plot(xl_bed, df_rmats_const_5ss[col_bed], fai, window)
    df_const_5ss = df_const_5ss.rename(columns={'coverage': 'const'})
    df_enhrest_3ss, df_coverage_enhrest_3ss, df_raw_enhrest_3ss = get_coverage_plot(xl_bed, df_rmats_enhrest_3ss[col_bed], fai, window)
    df_enhrest_3ss = df_enhrest_3ss.rename(columns={'coverage': 'enhanced_rest'})
    df_enhrest_5ss, df_coverage_enhrest_5ss, df_raw_enhrest_5ss = get_coverage_plot(xl_bed, df_rmats_enhrest_5ss[col_bed], fai, window)
    df_enhrest_5ss = df_enhrest_5ss.rename(columns={'coverage': 'enhanced_rest'})
    
    df_silrest_3ss, df_coverage_silrest_3ss,df_raw_silrest_3ss = get_coverage_plot(xl_bed, df_rmats_silrest_3ss[col_bed], fai, window)
    df_silrest_3ss = df_silrest_3ss.rename(columns={'coverage': 'silenced_rest'})
    df_silrest_5ss, df_coverage_silrest_5ss, df_raw_silrest_5ss = get_coverage_plot(xl_bed, df_rmats_silrest_5ss[col_bed], fai, window)
    df_silrest_5ss = df_silrest_5ss.rename(columns={'coverage': 'silenced_rest'})
    
    df_raw_enh_3ss['not_covered'] = len(df_rmats_enh_3ss) - df_raw_enh_3ss.coverage
    df_raw_enh_5ss['not_covered'] = len(df_rmats_enh_5ss) - df_raw_enh_5ss.coverage
    df_raw_sil_3ss['not_covered'] = len(df_rmats_sil_3ss) - df_raw_sil_3ss.coverage
    df_raw_sil_5ss['not_covered'] = len(df_rmats_sil_5ss) - df_raw_sil_5ss.coverage
    df_raw_ctrl_3ss['not_covered_control'] = len(df_rmats_ctrl_3ss) - df_raw_ctrl_3ss.coverage
    df_raw_ctrl_5ss['not_covered_control'] = len(df_rmats_ctrl_5ss) - df_raw_ctrl_5ss.coverage
    df_raw_enhrest_3ss['not_covered'] = len(df_rmats_enhrest_3ss) - df_raw_enhrest_3ss.coverage
    df_raw_enhrest_5ss['not_covered'] = len(df_rmats_enhrest_5ss) - df_raw_enhrest_5ss.coverage
    df_raw_silrest_3ss['not_covered'] = len(df_rmats_silrest_3ss) - df_raw_silrest_3ss.coverage
    df_raw_silrest_5ss['not_covered'] = len(df_rmats_silrest_5ss) - df_raw_silrest_5ss.coverage
    df_raw_const_3ss['not_covered'] = len(df_rmats_const_3ss) - df_raw_const_3ss.coverage
    df_raw_const_5ss['not_covered'] = len(df_rmats_const_5ss) - df_raw_const_5ss.coverage
    
    
    df_raw_ctrl_3ss = df_raw_ctrl_3ss.rename(columns={'coverage': 'control'})
    df_raw_ctrl_5ss = df_raw_ctrl_5ss.rename(columns={'coverage': 'control'})
    
    df_raw_enh_3ss = df_raw_enh_3ss.join(df_raw_ctrl_3ss)
    df_raw_enh_5ss = df_raw_enh_5ss.join(df_raw_ctrl_5ss)
    df_raw_sil_3ss = df_raw_sil_3ss.join(df_raw_ctrl_3ss)
    df_raw_sil_5ss = df_raw_sil_5ss.join(df_raw_ctrl_5ss)
    df_raw_enhrest_3ss = df_raw_enhrest_3ss.join(df_raw_ctrl_3ss)
    df_raw_enhrest_5ss = df_raw_enhrest_5ss.join(df_raw_ctrl_5ss)
    df_raw_silrest_3ss = df_raw_silrest_3ss.join(df_raw_ctrl_3ss)
    df_raw_silrest_5ss = df_raw_silrest_5ss.join(df_raw_ctrl_5ss)
    df_raw_const_3ss = df_raw_const_3ss.join(df_raw_ctrl_3ss)
    df_raw_const_5ss = df_raw_const_5ss.join(df_raw_ctrl_5ss)
       
    contingency_tables_enh_3ss = {}
    for index, row in df_raw_enh_3ss.iterrows():
        contingency_tables_enh_3ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_enh_5ss = {}
    for index, row in df_raw_enh_5ss.iterrows():
        contingency_tables_enh_5ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_sil_3ss = {}
    for index, row in df_raw_sil_3ss.iterrows():
        contingency_tables_sil_3ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_sil_5ss = {}
    for index, row in df_raw_sil_5ss.iterrows():
        contingency_tables_sil_5ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]      
    contingency_tables_enhrest_3ss = {}
    for index, row in df_raw_enhrest_3ss.iterrows():
        contingency_tables_enhrest_3ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_enhrest_5ss = {}
    for index, row in df_raw_enhrest_5ss.iterrows():
        contingency_tables_enhrest_5ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_silrest_3ss = {}
    for index, row in df_raw_silrest_3ss.iterrows():
        contingency_tables_silrest_3ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_silrest_5ss = {}
    for index, row in df_raw_silrest_5ss.iterrows():
        contingency_tables_silrest_5ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_const_3ss = {}
    for index, row in df_raw_const_3ss.iterrows():
        contingency_tables_const_3ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
    contingency_tables_const_5ss = {}
    for index, row in df_raw_const_5ss.iterrows():
        contingency_tables_const_5ss[index] = [[row['coverage'], row['not_covered']], 
                                              [row['control'], row['not_covered_control']]]
        
    fisher_result_enh_3ss = {}
    for k, v in contingency_tables_enh_3ss.items():
        fisher_result_enh_3ss[k] = stats.fisher_exact(v)
    fisher_result_enh_5ss = {}
    for k, v in contingency_tables_enh_5ss.items():
        fisher_result_enh_5ss[k] = stats.fisher_exact(v)
    fisher_result_sil_3ss = {}
    for k, v in contingency_tables_sil_3ss.items():
        fisher_result_sil_3ss[k] = stats.fisher_exact(v)
    fisher_result_sil_5ss = {}
    for k, v in contingency_tables_sil_5ss.items():
        fisher_result_sil_5ss[k] = stats.fisher_exact(v)
    fisher_result_enhrest_3ss = {}
    for k, v in contingency_tables_enhrest_3ss.items():
        fisher_result_enhrest_3ss[k] = stats.fisher_exact(v)
    fisher_result_enhrest_5ss = {}
    for k, v in contingency_tables_enhrest_5ss.items():
        fisher_result_enhrest_5ss[k] = stats.fisher_exact(v)
    fisher_result_silrest_3ss = {}
    for k, v in contingency_tables_silrest_3ss.items():
        fisher_result_silrest_3ss[k] = stats.fisher_exact(v)
    fisher_result_silrest_5ss = {}
    for k, v in contingency_tables_silrest_5ss.items():
        fisher_result_silrest_5ss[k] = stats.fisher_exact(v)
    fisher_result_const_3ss = {}
    for k, v in contingency_tables_const_3ss.items():
        fisher_result_const_3ss[k] = stats.fisher_exact(v)
    fisher_result_const_5ss = {}
    for k, v in contingency_tables_const_5ss.items():
        fisher_result_const_5ss[k] = stats.fisher_exact(v)
    p_values_fisher_result_enh_3ss = {}
    for k, v in fisher_result_enh_3ss.items():
        p_values_fisher_result_enh_3ss[k] = v[1]
    p_values_fisher_result_enh_5ss = {}
    for k, v in fisher_result_enh_5ss.items():
        p_values_fisher_result_enh_5ss[k] = v[1]
    p_values_fisher_result_sil_3ss = {}
    for k, v in fisher_result_sil_3ss.items():
        p_values_fisher_result_sil_3ss[k] = v[1]
    p_values_fisher_result_sil_5ss = {}
    for k, v in fisher_result_sil_5ss.items():
        p_values_fisher_result_sil_5ss[k] = v[1]
    p_values_fisher_result_enhrest_3ss = {}
    for k, v in fisher_result_enhrest_3ss.items():
        p_values_fisher_result_enhrest_3ss[k] = v[1]
    p_values_fisher_result_enhrest_5ss = {}
    for k, v in fisher_result_enhrest_5ss.items():
        p_values_fisher_result_enhrest_5ss[k] = v[1]
    p_values_fisher_result_silrest_3ss = {}
    for k, v in fisher_result_silrest_3ss.items():
        p_values_fisher_result_silrest_3ss[k] = v[1]
    p_values_fisher_result_silrest_5ss = {}
    for k, v in fisher_result_silrest_5ss.items():
        p_values_fisher_result_silrest_5ss[k] = v[1]
    p_values_fisher_result_const_3ss = {}
    for k, v in fisher_result_const_3ss.items():
        p_values_fisher_result_const_3ss[k] = v[1]
    p_values_fisher_result_const_5ss = {}
    for k, v in fisher_result_const_5ss.items():
        p_values_fisher_result_const_5ss[k] = v[1]
        
    df_fisher_3ss_enh = pd.DataFrame.from_dict(
        p_values_fisher_result_enh_3ss, orient='index', columns=['pval_enh'])
    df_fisher_3ss_sil = pd.DataFrame.from_dict(
        p_values_fisher_result_sil_3ss, orient='index', columns=['pval_sil'])
    df_fisher_5ss_enh = pd.DataFrame.from_dict(
        p_values_fisher_result_enh_5ss, orient='index', columns=['pval_enh'])
    df_fisher_5ss_sil = pd.DataFrame.from_dict(
        p_values_fisher_result_sil_5ss, orient='index', columns=['pval_sil'])
    df_fisher_3ss_enhrest = pd.DataFrame.from_dict(
        p_values_fisher_result_enhrest_3ss, orient='index', columns=['pval_enhrest'])
    df_fisher_3ss_silrest = pd.DataFrame.from_dict(
        p_values_fisher_result_silrest_3ss, orient='index', columns=['pval_silrest'])
    df_fisher_3ss_const = pd.DataFrame.from_dict(
        p_values_fisher_result_const_3ss, orient='index', columns=['pval_const'])
    df_fisher_5ss_enhrest = pd.DataFrame.from_dict(
        p_values_fisher_result_enhrest_5ss, orient='index', columns=['pval_enhrest'])
    df_fisher_5ss_silrest = pd.DataFrame.from_dict(
        p_values_fisher_result_silrest_5ss, orient='index', columns=['pval_silrest'])
    df_fisher_5ss_const = pd.DataFrame.from_dict(
        p_values_fisher_result_const_5ss, orient='index', columns=['pval_const'])
    
    df_fisher_3ss = df_fisher_3ss_enh.join([
        df_fisher_3ss_sil, df_enh_3ss, df_sil_3ss, df_ctrl_3ss,
        df_fisher_3ss_enhrest, df_fisher_3ss_silrest, df_fisher_3ss_const,
        df_enhrest_3ss, df_silrest_3ss, df_const_3ss])
    df_fisher_5ss = df_fisher_5ss_enh.join([
        df_fisher_5ss_sil, df_enh_5ss, df_sil_5ss, df_ctrl_5ss,
        df_fisher_5ss_enhrest, df_fisher_5ss_silrest, df_fisher_5ss_const,
        df_enhrest_5ss, df_silrest_5ss, df_const_5ss])
    
    df_fisher_3ss['enh_fold_change'] = df_fisher_3ss.enhanced / df_fisher_3ss.control
    df_fisher_3ss['sil_fold_change'] = df_fisher_3ss.silenced / df_fisher_3ss.control
    df_fisher_5ss['enh_fold_change'] = df_fisher_5ss.enhanced / df_fisher_5ss.control
    df_fisher_5ss['sil_fold_change'] = df_fisher_5ss.silenced / df_fisher_5ss.control    
    df_fisher_3ss['enhrest_fold_change'] = df_fisher_3ss.enhanced_rest / df_fisher_3ss.control
    df_fisher_3ss['silrest_fold_change'] = df_fisher_3ss.silenced_rest / df_fisher_3ss.control
    df_fisher_3ss['const_fold_change'] = df_fisher_3ss.const / df_fisher_3ss.control
    df_fisher_5ss['enhrest_fold_change'] = df_fisher_5ss.enhanced_rest / df_fisher_5ss.control
    df_fisher_5ss['silrest_fold_change'] = df_fisher_5ss.silenced_rest / df_fisher_5ss.control
    df_fisher_5ss['const_fold_change'] = df_fisher_5ss.const / df_fisher_5ss.control
    
    df_fisher_3ss['-log(pval)_enh'] = -np.log(df_fisher_3ss.pval_enh)
    df_fisher_3ss['-log(pval)_sil'] = -np.log(df_fisher_3ss.pval_sil)
    df_fisher_5ss['-log(pval)_enh'] = -np.log(df_fisher_5ss.pval_enh)
    df_fisher_5ss['-log(pval)_sil'] = -np.log(df_fisher_5ss.pval_sil)    
    df_fisher_3ss['-log(pval)_enhrest'] = -np.log(df_fisher_3ss.pval_enhrest)
    df_fisher_3ss['-log(pval)_silrest'] = -np.log(df_fisher_3ss.pval_silrest)
    df_fisher_3ss['-log(pval)_const'] = -np.log(df_fisher_3ss.pval_const)
    df_fisher_5ss['-log(pval)_enhrest'] = -np.log(df_fisher_5ss.pval_enhrest)
    df_fisher_5ss['-log(pval)_silrest'] = -np.log(df_fisher_5ss.pval_silrest)
    df_fisher_5ss['-log(pval)_const'] = -np.log(df_fisher_5ss.pval_const)
    
    df_fisher_3ss['-log(pval)_enh_forplot'] = np.nan
    df_fisher_3ss['-log(pval)_sil_forplot'] = np.nan
    df_fisher_3ss['-log(pval)_enhrest_forplot'] = np.nan
    df_fisher_3ss['-log(pval)_silrest_forplot'] = np.nan
    df_fisher_3ss['-log(pval)_const_forplot'] = np.nan
    for index, row in df_fisher_3ss.iterrows():
        if row['enh_fold_change'] > 1:
            df_fisher_3ss.loc[index, '-log(pval)_enh_forplot'] = row['-log(pval)_enh']
        else:
            df_fisher_3ss.loc[index, '-log(pval)_enh_forplot'] = -row['-log(pval)_enh']
        if row['sil_fold_change'] > 1:
            df_fisher_3ss.loc[index, '-log(pval)_sil_forplot'] = row['-log(pval)_sil']
        else:
            df_fisher_3ss.loc[index, '-log(pval)_sil_forplot'] = -row['-log(pval)_sil']
        if row['enhrest_fold_change'] > 1:
            df_fisher_3ss.loc[index, '-log(pval)_enhrest_forplot'] = row['-log(pval)_enhrest']
        else:
            df_fisher_3ss.loc[index, '-log(pval)_enhrest_forplot'] = -row['-log(pval)_enhrest']
        if row['silrest_fold_change'] > 1:
            df_fisher_3ss.loc[index, '-log(pval)_silrest_forplot'] = row['-log(pval)_silrest']
        else:
            df_fisher_3ss.loc[index, '-log(pval)_silrest_forplot'] = -row['-log(pval)_silrest']
        if row['const_fold_change'] > 1:
            df_fisher_3ss.loc[index, '-log(pval)_const_forplot'] = row['-log(pval)_const']
        else:
            df_fisher_3ss.loc[index, '-log(pval)_const_forplot'] = -row['-log(pval)_const']
    
    df_fisher_5ss['-log(pval)_enh_forplot'] = np.nan
    df_fisher_5ss['-log(pval)_sil_forplot'] = np.nan
    df_fisher_5ss['-log(pval)_enhrest_forplot'] = np.nan
    df_fisher_5ss['-log(pval)_silrest_forplot'] = np.nan
    df_fisher_5ss['-log(pval)_const_forplot'] = np.nan
    for index, row in df_fisher_5ss.iterrows():
        if row['enh_fold_change'] > 1:
            df_fisher_5ss.loc[index, '-log(pval)_enh_forplot'] = row['-log(pval)_enh']
        else:
            df_fisher_5ss.loc[index, '-log(pval)_enh_forplot'] = -row['-log(pval)_enh']
        if row['sil_fold_change'] > 1:
            df_fisher_5ss.loc[index, '-log(pval)_sil_forplot'] = row['-log(pval)_sil']
        else:
            df_fisher_5ss.loc[index, '-log(pval)_sil_forplot'] = -row['-log(pval)_sil']
        if row['enhrest_fold_change'] > 1:
            df_fisher_5ss.loc[index, '-log(pval)_enhrest_forplot'] = row['-log(pval)_enhrest']
        else:
            df_fisher_5ss.loc[index, '-log(pval)_enhrest_forplot'] = -row['-log(pval)_enhrest']
        if row['silrest_fold_change'] > 1:
            df_fisher_5ss.loc[index, '-log(pval)_silrest_forplot'] = row['-log(pval)_silrest']
        else:
            df_fisher_5ss.loc[index, '-log(pval)_silrest_forplot'] = -row['-log(pval)_silrest']
        if row['const_fold_change'] > 1:
            df_fisher_5ss.loc[index, '-log(pval)_const_forplot'] = row['-log(pval)_const']
        else:
            df_fisher_5ss.loc[index, '-log(pval)_const_forplot'] = -row['-log(pval)_const']
    
    df_fisher_3ss['enh_smooth'] = df_fisher_3ss['-log(pval)_enh_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enh_max'] = df_fisher_3ss['-log(pval)_enh_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['enh_min'] = df_fisher_3ss['-log(pval)_enh_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['enh_smooth'] = df_fisher_5ss['-log(pval)_enh_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enh_max'] = df_fisher_5ss['-log(pval)_enh_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['enh_min'] = df_fisher_5ss['-log(pval)_enh_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['enhrest_smooth'] = df_fisher_3ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enhrest_max'] = df_fisher_3ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['enhrest_min'] = df_fisher_3ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['enhrest_smooth'] = df_fisher_5ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enhrest_max'] = df_fisher_5ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['enhrest_min'] = df_fisher_5ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['sil_smooth'] = df_fisher_3ss['-log(pval)_sil_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['sil_max'] = df_fisher_3ss['-log(pval)_sil_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['sil_min'] = df_fisher_3ss['-log(pval)_sil_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['sil_smooth'] = df_fisher_5ss['-log(pval)_sil_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['sil_max'] = df_fisher_5ss['-log(pval)_sil_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['sil_min'] = df_fisher_5ss['-log(pval)_sil_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['silrest_smooth'] = df_fisher_3ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['silrest_max'] = df_fisher_3ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['silrest_min'] = df_fisher_3ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['silrest_smooth'] = df_fisher_5ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['silrest_max'] = df_fisher_5ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['silrest_min'] = df_fisher_5ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['enh_smooth_real'] = df_fisher_3ss['-log(pval)_enh'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enh_max_real'] = df_fisher_3ss['-log(pval)_enh'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['enh_min_real'] = df_fisher_3ss['-log(pval)_enh'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['enh_smooth_real'] = df_fisher_5ss['-log(pval)_enh'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enh_max_real'] = df_fisher_5ss['-log(pval)_enh'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['enh_min_real'] = df_fisher_5ss['-log(pval)_enh'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    
    df_fisher_3ss['enhrest_smooth_real'] = df_fisher_3ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enhrest_max_real'] = df_fisher_3ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['enhrest_min_real'] = df_fisher_3ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['enhrest_smooth_real'] = df_fisher_5ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enhrest_max_real'] = df_fisher_5ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['enhrest_min_real'] = df_fisher_5ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['sil_smooth_real'] = df_fisher_3ss['-log(pval)_sil'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['sil_max_real'] = df_fisher_3ss['-log(pval)_sil'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['sil_min_real'] = df_fisher_3ss['-log(pval)_sil'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['sil_smooth_real'] = df_fisher_5ss['-log(pval)_sil'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['sil_max_real'] = df_fisher_5ss['-log(pval)_sil'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['sil_min_real'] = df_fisher_5ss['-log(pval)_sil'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)

    df_fisher_3ss['enhrest_smooth'] = df_fisher_3ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enhrest_max'] = df_fisher_3ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['enhrest_min'] = df_fisher_3ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['enhrest_smooth'] = df_fisher_5ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enhrest_max'] = df_fisher_5ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['enhrest_min'] = df_fisher_5ss['-log(pval)_enhrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['silrest_smooth'] = df_fisher_3ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['silrest_max'] = df_fisher_3ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['silrest_min'] = df_fisher_3ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['silrest_smooth'] = df_fisher_5ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['silrest_max'] = df_fisher_5ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['silrest_min'] = df_fisher_5ss['-log(pval)_silrest_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['const_smooth'] = df_fisher_3ss['-log(pval)_const_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['const_max'] = df_fisher_3ss['-log(pval)_const_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['const_min'] = df_fisher_3ss['-log(pval)_const_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['const_smooth'] = df_fisher_5ss['-log(pval)_const_forplot'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['const_max'] = df_fisher_5ss['-log(pval)_const_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['const_min'] = df_fisher_5ss['-log(pval)_const_forplot'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['enhrest_smooth_real'] = df_fisher_3ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enhrest_max_real'] = df_fisher_3ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['enhrest_min_real'] = df_fisher_3ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['enhrest_smooth_real'] = df_fisher_5ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enhrest_max_real'] = df_fisher_5ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['enhrest_min_real'] = df_fisher_5ss['-log(pval)_enhrest'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['silrest_smooth_real'] = df_fisher_3ss['-log(pval)_silrest'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['silrest_max_real'] = df_fisher_3ss['-log(pval)_silrest'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['silrest_min_real'] = df_fisher_3ss['-log(pval)_silrest'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['silrest_smooth_real'] = df_fisher_5ss['-log(pval)_silrest'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['silrest_max_real'] = df_fisher_5ss['-log(pval)_silrest'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['silrest_min_real'] = df_fisher_5ss['-log(pval)_silrest'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    
    df_fisher_3ss['const_smooth_real'] = df_fisher_3ss['-log(pval)_const'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['const_max_real'] = df_fisher_3ss['-log(pval)_const'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_3ss['const_min_real'] = df_fisher_3ss['-log(pval)_const'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)
    df_fisher_5ss['const_smooth_real'] = df_fisher_5ss['-log(pval)_const'].rolling(
        smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['const_max_real'] = df_fisher_5ss['-log(pval)_const'].rolling(
        smoothing, center=True).apply(lambda x: np.max(x), raw=False)
    df_fisher_5ss['const_min_real'] = df_fisher_5ss['-log(pval)_const'].rolling(
        smoothing, center=True).apply(lambda x: np.min(x), raw=False)

    
    enh_3ss_no_xl = get_exon_dist(len(df_rmats_enh_3ss), df_coverage_enh_3ss, window)
    enh_5ss_no_xl = get_exon_dist(len(df_rmats_enh_5ss), df_coverage_enh_5ss, window)
    enh_no_xl = {k: enh_3ss_no_xl[k] + enh_5ss_no_xl[k] for k in enh_3ss_no_xl.keys()}
    sil_3ss_no_xl = get_exon_dist(len(df_rmats_sil_3ss), df_coverage_sil_3ss, window)
    sil_5ss_no_xl = get_exon_dist(len(df_rmats_sil_3ss), df_coverage_sil_5ss, window)
    sil_no_xl = {k: sil_3ss_no_xl[k] + sil_5ss_no_xl[k] for k in sil_3ss_no_xl.keys()}
    
    enhrest_3ss_no_xl = get_exon_dist(len(df_rmats_enhrest_3ss), df_coverage_enhrest_3ss, window)
    enhrest_5ss_no_xl = get_exon_dist(len(df_rmats_enhrest_5ss), df_coverage_enhrest_5ss, window)
    enhrest_no_xl = {k: enhrest_3ss_no_xl[k] + enhrest_5ss_no_xl[k] for k in enhrest_3ss_no_xl.keys()}
    silrest_3ss_no_xl = get_exon_dist(len(df_rmats_silrest_3ss), df_coverage_silrest_3ss, window)
    silrest_5ss_no_xl = get_exon_dist(len(df_rmats_silrest_3ss), df_coverage_silrest_5ss, window)
    silrest_no_xl = {k: silrest_3ss_no_xl[k] + silrest_5ss_no_xl[k] for k in silrest_3ss_no_xl.keys()}
    
    const_3ss_no_xl = get_exon_dist(len(df_rmats_const_3ss), df_coverage_const_3ss, window)
    const_5ss_no_xl = get_exon_dist(len(df_rmats_const_3ss), df_coverage_const_5ss, window)
    const_no_xl = {k: const_3ss_no_xl[k] + const_5ss_no_xl[k] for k in const_3ss_no_xl.keys()}
                               
    ctrl_3ss_no_xl = get_exon_dist(len(df_rmats_ctrl_3ss), df_coverage_ctrl_3ss, window)
    ctrl_5ss_no_xl = get_exon_dist(len(df_rmats_ctrl_3ss), df_coverage_ctrl_5ss, window)
    ctrl_no_xl = {k: ctrl_3ss_no_xl[k] + ctrl_5ss_no_xl[k] for k in ctrl_3ss_no_xl.keys()}
    
    if de_source == 'rmats':
        col_upstream = ['chr', 'upstream_es', 'upstream_ee', 'name', 'score', 'strand']
        upstream_rename = {'upstream_es': 'start', 'upstream_ee': 'end'}
        df_enh_3ss_upstream, _, _ = get_coverage_plot(xl_bed, df_rmats_enh_3ss[col_upstream].rename(
            columns=upstream_rename), fai, window)
        df_enh_3ss_upstream = df_enh_3ss_upstream.rename(columns={'coverage': 'enhanced'})
        df_enh_5ss_upstream, _, _ = get_coverage_plot(xl_bed, df_rmats_enh_5ss[col_upstream].rename(
            columns=upstream_rename), fai, window)
        df_enh_5ss_upstream = df_enh_5ss_upstream.rename(columns={'coverage': 'enhanced'})
        df_sil_3ss_upstream, _ , _ = get_coverage_plot(xl_bed, df_rmats_sil_3ss[col_upstream].rename(
            columns=upstream_rename), fai, window)
        df_sil_3ss_upstream = df_sil_3ss_upstream.rename(columns={'coverage': 'silenced'})
        df_sil_5ss_upstream, _, _ = get_coverage_plot(xl_bed, df_rmats_sil_5ss[col_upstream].rename(
            columns=upstream_rename), fai, window)
        df_sil_5ss_upstream = df_sil_5ss_upstream.rename(columns={'coverage': 'silenced'})                                                           
        df_ctrl_3ss_upstream, _, _ = get_coverage_plot(xl_bed, df_rmats_ctrl_3ss[col_upstream].rename(
            columns=upstream_rename), fai, window)
        df_ctrl_3ss_upstream = df_ctrl_3ss_upstream.rename(columns={'coverage': 'control'})
        df_ctrl_5ss_upstream, _, _ = get_coverage_plot(xl_bed, df_rmats_ctrl_5ss[col_upstream].rename(
            columns=upstream_rename), fai, window)
        df_ctrl_5ss_upstream = df_ctrl_5ss_upstream.rename(columns={'coverage': 'control'})

        col_downstream = ['chr', 'downstream_es', 'downstream_ee', 'name', 'score', 'strand']
        downstream_rename = {'downstream_es': 'start', 'downstream_ee': 'end'}
        df_enh_3ss_downstream, _, _ = get_coverage_plot(xl_bed, df_rmats_enh_3ss[col_downstream].rename(
            columns=downstream_rename), fai, window)
        df_enh_3ss_downstream = df_enh_3ss_downstream.rename(columns={'coverage': 'enhanced'})
        df_enh_5ss_downstream, _, _ = get_coverage_plot(xl_bed, df_rmats_enh_5ss[col_downstream].rename(
            columns=downstream_rename), fai, window)
        df_enh_5ss_downstream = df_enh_5ss_downstream.rename(columns={'coverage': 'enhanced'})
        df_sil_3ss_downstream, _, _ = get_coverage_plot(xl_bed, df_rmats_sil_3ss[col_downstream].rename(
            columns=downstream_rename), fai, window)
        df_sil_3ss_downstream = df_sil_3ss_downstream.rename(columns={'coverage': 'silenced'})
        df_sil_5ss_downstream, _, _ = get_coverage_plot(xl_bed, df_rmats_sil_5ss[col_downstream].rename(
            columns=downstream_rename), fai, window)
        df_sil_5ss_downstream = df_sil_5ss_downstream.rename(columns={'coverage': 'silenced'})                                                          
        df_ctrl_3ss_downstream, _, _ = get_coverage_plot(xl_bed, df_rmats_ctrl_3ss[col_downstream].rename(
            columns=downstream_rename), fai, window)
        df_ctrl_3ss_downstream = df_ctrl_3ss_downstream.rename(columns={'coverage': 'control'})
        df_ctrl_5ss_downstream, _, _ = get_coverage_plot(xl_bed, df_rmats_ctrl_5ss[col_downstream].rename(
            columns=downstream_rename), fai, window)
        df_ctrl_5ss_downstream = df_ctrl_5ss_downstream.rename(columns={'coverage': 'control'})

    df_final_3ss = df_enh_3ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2).merge(
        df_sil_3ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_ctrl_3ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_enhrest_3ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_silrest_3ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_const_3ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True)
    
    df_final_5ss = df_enh_5ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2).merge(
        df_sil_5ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_ctrl_5ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_enhrest_5ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_silrest_5ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
        df_const_5ss.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True)

    df_final_3ss.to_csv(f'{output_dir}/{name}_final_3ss.tsv', sep='\t', index=None)
    df_final_5ss.to_csv(f'{output_dir}/{name}_final_5ss.tsv', sep='\t', index=None)


    df_rmats_enh_3ss_temp = df_rmats_enh_3ss.reset_index()
    rmats_enh_3ss = pbt.BedTool.from_dataframe(df_rmats_enh_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enh_3ss = rmats_enh_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enh_3ss = rmats_enh_3ss.to_dataframe(header=None)
    df_enh_3ss = df_enh_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enh_3ss_grouped10 = df_enh_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enh_3ss_grouped_12 = df_enh_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enh_3ss = df_enh_3ss_grouped10.join(df_enh_3ss_grouped_12).reset_index()
    df_enh_3ss = df_enh_3ss.set_index(3)
    df_rmats_enh_3ss_temp = pd.merge(df_rmats_enh_3ss_temp, df_enh_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enh_3ss_temp['index']
    df_rmats_enh_3ss_temp = df_rmats_enh_3ss_temp.set_index(14)
    df_rmats_enh_3ss_temp = df_rmats_enh_3ss.reset_index()
    rmats_enh_3ss = pbt.BedTool.from_dataframe(df_rmats_enh_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enh_3ss = rmats_enh_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enh_3ss = rmats_enh_3ss.to_dataframe(header=None)
    df_enh_3ss = df_enh_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enh_3ss_grouped10 = df_enh_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enh_3ss_grouped_12 = df_enh_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enh_3ss = df_enh_3ss_grouped10.join(df_enh_3ss_grouped_12).reset_index()
    df_enh_3ss = df_enh_3ss.set_index(3)
    df_rmats_enh_3ss_temp = pd.merge(df_rmats_enh_3ss_temp, df_enh_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enh_3ss_temp['index']
    df_rmats_enh_3ss_temp = df_rmats_enh_3ss_temp.set_index(14)
    df_rmats_enh_3ss_temp = df_rmats_enh_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss.reset_index()
    rmats_enh_5ss = pbt.BedTool.from_dataframe(df_rmats_enh_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enh_5ss = rmats_enh_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enh_5ss = rmats_enh_5ss.to_dataframe(header=None)
    df_enh_5ss = df_enh_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enh_5ss_grouped10 = df_enh_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enh_5ss_grouped_12 = df_enh_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enh_5ss = df_enh_5ss_grouped10.join(df_enh_5ss_grouped_12).reset_index()
    df_enh_5ss = df_enh_5ss.set_index(3)
    df_rmats_enh_5ss_temp = pd.merge(df_rmats_enh_5ss_temp, df_enh_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enh_5ss_temp['index']
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss_temp.set_index(14)
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_enh = df_rmats_enh_5ss_temp.join(df_rmats_enh_3ss_temp, lsuffix='5ss', rsuffix='3ss')


    df_rmats_enh2 = df_rmats_enh[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
       'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_enh2 = df_rmats_enh2.reset_index()

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & (2*window > (df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss'])) & \
                      ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) > window),
                      'start'] = df_rmats_enh2['start5ss'] - window
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & (2*window > (df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss'])) & \
                      ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) > window),
                      'end'] = df_rmats_enh2['start3ss'] + window

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & (2*window > (df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss'])) & \
                      ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) > window),
                      'start'] = df_rmats_enh2['start3ss'] - window
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & (2*window > (df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss'])) & \
                      ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) > window),
                      'end'] = df_rmats_enh2['start5ss'] + window

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) < window),
                      'start'] = df_rmats_enh2['start3ss']
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) < window),
                      'end'] = df_rmats_enh2['start5ss']

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) < window),
                      'start'] = df_rmats_enh2['start5ss']
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) < window),
                      'end'] = df_rmats_enh2['start3ss']

    df_rmats_enh2['len'] = df_rmats_enh2['end'] - df_rmats_enh2['start']

    df_rmats_enh2['diff'] = abs(df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss'])

    df_rmats_enh2 = df_rmats_enh2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_enh2 = df_rmats_enh2[~(df_rmats_enh2['start'].isna() | df_rmats_enh2['end'].isna())]

    df_rmats_enh2['start'] = df_rmats_enh2['start'].astype(int)
    df_rmats_enh2['end'] = df_rmats_enh2['end'].astype(int)

    rmats_enh2 = pbt.BedTool.from_dataframe(df_rmats_enh2)
    rmats_enh2 = rmats_enh2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_enh2 = rmats_enh2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_enh2['base_coverage_5ss']
    del df_rmats_enh2['chrom2']
    del df_rmats_enh2['start2']
    del df_rmats_enh2['end2']
    del df_rmats_enh2['useless']
    del df_rmats_enh2['strand2']

    df_rmats_enh2_grouped10 = df_rmats_enh2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_enh2_grouped12 = df_rmats_enh2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_enh2 = df_rmats_enh2_grouped10.join(df_rmats_enh2_grouped12).reset_index()
    df_rmats_enh2 = df_rmats_enh2.set_index('index')

    df_rmats_enh_overlap = df_rmats_enh.join(df_rmats_enh2)

    df_rmats_enh_overlap = df_rmats_enh_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_enh_overlap.loc[df_rmats_enh_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_enh_overlap.loc[df_rmats_enh_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_enh_overlap.loc[df_rmats_enh_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_enh_overlap['cDNA_coverage_5ss'] = df_rmats_enh_overlap['cDNA_coverage_5ss'] - df_rmats_enh_overlap['overlap_score']
    df_rmats_enh_overlap['cDNA_coverage_3ss'] = df_rmats_enh_overlap['cDNA_coverage_3ss'] - df_rmats_enh_overlap['overlap_score']
    df_rmats_enh_overlap['base_coverage_5ss'] = df_rmats_enh_overlap['base_coverage_5ss'] - df_rmats_enh_overlap['overlap_base']
    df_rmats_enh_overlap['base_coverage_3ss'] = df_rmats_enh_overlap['base_coverage_3ss'] - df_rmats_enh_overlap['overlap_base']

    df_rmats_enh = df_rmats_enh_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_sil_3ss_temp = df_rmats_sil_3ss.reset_index()
    rmats_sil_3ss = pbt.BedTool.from_dataframe(df_rmats_sil_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_sil_3ss = rmats_sil_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_sil_3ss = rmats_sil_3ss.to_dataframe(header=None)
    df_sil_3ss = df_sil_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_sil_3ss_grouped10 = df_sil_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_sil_3ss_grouped_12 = df_sil_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_sil_3ss = df_sil_3ss_grouped10.join(df_sil_3ss_grouped_12).reset_index()
    df_sil_3ss = df_sil_3ss.set_index(3)
    df_rmats_sil_3ss_temp = pd.merge(df_rmats_sil_3ss_temp, df_sil_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_sil_3ss_temp['index']
    df_rmats_sil_3ss_temp = df_rmats_sil_3ss_temp.set_index(14)
    df_rmats_sil_3ss_temp = df_rmats_sil_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_sil_5ss_temp = df_rmats_sil_5ss.reset_index()
    rmats_sil_5ss = pbt.BedTool.from_dataframe(df_rmats_sil_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_sil_5ss = rmats_sil_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_sil_5ss = rmats_sil_5ss.to_dataframe(header=None)
    df_sil_5ss = df_sil_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_sil_5ss_grouped10 = df_sil_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_sil_5ss_grouped_12 = df_sil_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_sil_5ss = df_sil_5ss_grouped10.join(df_sil_5ss_grouped_12).reset_index()
    df_sil_5ss = df_sil_5ss.set_index(3)
    df_rmats_sil_5ss_temp = pd.merge(df_rmats_sil_5ss_temp, df_sil_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_sil_5ss_temp['index']
    df_rmats_sil_5ss_temp = df_rmats_sil_5ss_temp.set_index(14)
    df_rmats_sil_5ss_temp = df_rmats_sil_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_sil = df_rmats_sil_5ss_temp.join(df_rmats_sil_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_sil2 = df_rmats_sil[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
       'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_sil2 = df_rmats_sil2.reset_index()

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & (2*window > (df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss'])) & \
                      ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) > window),
                      'start'] = df_rmats_sil2['start5ss'] - window
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & (2*window > (df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss'])) & \
                      ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) > window),
                      'end'] = df_rmats_sil2['start3ss'] + window

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & (2*window > (df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss'])) & \
                      ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) > window),
                      'start'] = df_rmats_sil2['start3ss'] - window
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & (2*window > (df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss'])) & \
                      ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) > window),
                      'end'] = df_rmats_sil2['start5ss'] + window

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) < window),
                      'start'] = df_rmats_sil2['start3ss']
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) < window),
                      'end'] = df_rmats_sil2['start5ss']

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) < window),
                      'start'] = df_rmats_sil2['start5ss']
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) < window),
                      'end'] = df_rmats_sil2['start3ss']

    df_rmats_sil2['len'] = df_rmats_sil2['end'] - df_rmats_sil2['start']

    df_rmats_sil2['diff'] = abs(df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss'])

    df_rmats_sil2 = df_rmats_sil2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_sil2 = df_rmats_sil2[~(df_rmats_sil2['start'].isna() | df_rmats_sil2['end'].isna())]

    df_rmats_sil2['start'] = df_rmats_sil2['start'].astype(int)
    df_rmats_sil2['end'] = df_rmats_sil2['end'].astype(int)

    rmats_sil2 = pbt.BedTool.from_dataframe(df_rmats_sil2)
    rmats_sil2 = rmats_sil2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_sil2 = rmats_sil2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_sil2['base_coverage_5ss']
    del df_rmats_sil2['chrom2']
    del df_rmats_sil2['start2']
    del df_rmats_sil2['end2']
    del df_rmats_sil2['useless']
    del df_rmats_sil2['strand2']

    df_rmats_sil2_grouped10 = df_rmats_sil2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_sil2_grouped12 = df_rmats_sil2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_sil2 = df_rmats_sil2_grouped10.join(df_rmats_sil2_grouped12).reset_index()
    df_rmats_sil2 = df_rmats_sil2.set_index('index')

    df_rmats_sil_overlap = df_rmats_sil.join(df_rmats_sil2)

    df_rmats_sil_overlap = df_rmats_sil_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_sil_overlap.loc[df_rmats_sil_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_sil_overlap.loc[df_rmats_sil_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_sil_overlap.loc[df_rmats_sil_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_sil_overlap['cDNA_coverage_5ss'] = df_rmats_sil_overlap['cDNA_coverage_5ss'] - df_rmats_sil_overlap['overlap_score']
    df_rmats_sil_overlap['cDNA_coverage_3ss'] = df_rmats_sil_overlap['cDNA_coverage_3ss'] - df_rmats_sil_overlap['overlap_score']
    df_rmats_sil_overlap['base_coverage_5ss'] = df_rmats_sil_overlap['base_coverage_5ss'] - df_rmats_sil_overlap['overlap_base']
    df_rmats_sil_overlap['base_coverage_3ss'] = df_rmats_sil_overlap['base_coverage_3ss'] - df_rmats_sil_overlap['overlap_base']

    df_rmats_sil = df_rmats_sil_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss.reset_index()
    rmats_enhrest_3ss = pbt.BedTool.from_dataframe(df_rmats_enhrest_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enhrest_3ss = rmats_enhrest_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enhrest_3ss = rmats_enhrest_3ss.to_dataframe(header=None)
    df_enhrest_3ss = df_enhrest_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enhrest_3ss_grouped10 = df_enhrest_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enhrest_3ss_grouped_12 = df_enhrest_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enhrest_3ss = df_enhrest_3ss_grouped10.join(df_enhrest_3ss_grouped_12).reset_index()
    df_enhrest_3ss = df_enhrest_3ss.set_index(3)
    df_rmats_enhrest_3ss_temp = pd.merge(df_rmats_enhrest_3ss_temp, df_enhrest_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enhrest_3ss_temp['index']
    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss_temp.set_index(14)
    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss.reset_index()
    rmats_enhrest_5ss = pbt.BedTool.from_dataframe(df_rmats_enhrest_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enhrest_5ss = rmats_enhrest_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enhrest_5ss = rmats_enhrest_5ss.to_dataframe(header=None)
    df_enhrest_5ss = df_enhrest_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enhrest_5ss_grouped10 = df_enhrest_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enhrest_5ss_grouped_12 = df_enhrest_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enhrest_5ss = df_enhrest_5ss_grouped10.join(df_enhrest_5ss_grouped_12).reset_index()
    df_enhrest_5ss = df_enhrest_5ss.set_index(3)
    df_rmats_enhrest_5ss_temp = pd.merge(df_rmats_enhrest_5ss_temp, df_enhrest_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enhrest_5ss_temp['index']
    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss_temp.set_index(14)
    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_enhrest = df_rmats_enhrest_5ss_temp.join(df_rmats_enhrest_3ss_temp, lsuffix='5ss', rsuffix='3ss')

    df_rmats_enhrest2 = df_rmats_enhrest[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_enhrest2 = df_rmats_enhrest2.reset_index()

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & (2*window > (df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss'])) & \
                    ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) > window),
                    'start'] = df_rmats_enhrest2['start5ss'] - window
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & (2*window > (df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss'])) & \
                    ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) > window),
                    'end'] = df_rmats_enhrest2['start3ss'] + window

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & (2*window > (df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss'])) & \
                    ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) > window),
                    'start'] = df_rmats_enhrest2['start3ss'] - window
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & (2*window > (df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss'])) & \
                    ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) > window),
                    'end'] = df_rmats_enhrest2['start5ss'] + window

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) < window),
                    'start'] = df_rmats_enhrest2['start3ss']
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) < window),
                    'end'] = df_rmats_enhrest2['start5ss']

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) < window),
                    'start'] = df_rmats_enhrest2['start5ss']
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) < window),
                    'end'] = df_rmats_enhrest2['start3ss']

    df_rmats_enhrest2['len'] = df_rmats_enhrest2['end'] - df_rmats_enhrest2['start']

    df_rmats_enhrest2['diff'] = abs(df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss'])

    df_rmats_enhrest2 = df_rmats_enhrest2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_enhrest2 = df_rmats_enhrest2[~(df_rmats_enhrest2['start'].isna() | df_rmats_enhrest2['end'].isna())]

    df_rmats_enhrest2['start'] = df_rmats_enhrest2['start'].astype(int)
    df_rmats_enhrest2['end'] = df_rmats_enhrest2['end'].astype(int)

    rmats_enhrest2 = pbt.BedTool.from_dataframe(df_rmats_enhrest2)
    rmats_enhrest2 = rmats_enhrest2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_enhrest2 = rmats_enhrest2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_enhrest2['base_coverage_5ss']
    del df_rmats_enhrest2['chrom2']
    del df_rmats_enhrest2['start2']
    del df_rmats_enhrest2['end2']
    del df_rmats_enhrest2['useless']
    del df_rmats_enhrest2['strand2']

    df_rmats_enhrest2_grouped10 = df_rmats_enhrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_enhrest2_grouped12 = df_rmats_enhrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_enhrest2 = df_rmats_enhrest2_grouped10.join(df_rmats_enhrest2_grouped12).reset_index()
    df_rmats_enhrest2 = df_rmats_enhrest2.set_index('index')

    df_rmats_enhrest_overlap = df_rmats_enhrest.join(df_rmats_enhrest2)

    df_rmats_enhrest_overlap = df_rmats_enhrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_enhrest_overlap.loc[df_rmats_enhrest_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_enhrest_overlap.loc[df_rmats_enhrest_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_enhrest_overlap.loc[df_rmats_enhrest_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_enhrest_overlap['cDNA_coverage_5ss'] = df_rmats_enhrest_overlap['cDNA_coverage_5ss'] - df_rmats_enhrest_overlap['overlap_score']
    df_rmats_enhrest_overlap['cDNA_coverage_3ss'] = df_rmats_enhrest_overlap['cDNA_coverage_3ss'] - df_rmats_enhrest_overlap['overlap_score']
    df_rmats_enhrest_overlap['base_coverage_5ss'] = df_rmats_enhrest_overlap['base_coverage_5ss'] - df_rmats_enhrest_overlap['overlap_base']
    df_rmats_enhrest_overlap['base_coverage_3ss'] = df_rmats_enhrest_overlap['base_coverage_3ss'] - df_rmats_enhrest_overlap['overlap_base']

    df_rmats_enhrest = df_rmats_enhrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    
    
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss.reset_index()
    rmats_silrest_3ss = pbt.BedTool.from_dataframe(df_rmats_silrest_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_silrest_3ss = rmats_silrest_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_silrest_3ss = rmats_silrest_3ss.to_dataframe(header=None)
    df_silrest_3ss = df_silrest_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_silrest_3ss_grouped10 = df_silrest_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_silrest_3ss_grouped_12 = df_silrest_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_silrest_3ss = df_silrest_3ss_grouped10.join(df_silrest_3ss_grouped_12).reset_index()
    df_silrest_3ss = df_silrest_3ss.set_index(3)
    df_rmats_silrest_3ss_temp = pd.merge(df_rmats_silrest_3ss_temp, df_silrest_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_silrest_3ss_temp['index']
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss_temp.set_index(14)
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss.reset_index()
    rmats_silrest_5ss = pbt.BedTool.from_dataframe(df_rmats_silrest_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_silrest_5ss = rmats_silrest_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_silrest_5ss = rmats_silrest_5ss.to_dataframe(header=None)
    df_silrest_5ss = df_silrest_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_silrest_5ss_grouped10 = df_silrest_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_silrest_5ss_grouped_12 = df_silrest_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_silrest_5ss = df_silrest_5ss_grouped10.join(df_silrest_5ss_grouped_12).reset_index()
    df_silrest_5ss = df_silrest_5ss.set_index(3)
    df_rmats_silrest_5ss_temp = pd.merge(df_rmats_silrest_5ss_temp, df_silrest_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_silrest_5ss_temp['index']
    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss_temp.set_index(14)
    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_silrest = df_rmats_silrest_5ss_temp.join(df_rmats_silrest_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_silrest2 = df_rmats_silrest[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_silrest2 = df_rmats_silrest2.reset_index()

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & (2*window > (df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss'])) & \
                    ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) > window),
                    'start'] = df_rmats_silrest2['start5ss'] - window
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & (2*window > (df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss'])) & \
                    ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) > window),
                    'end'] = df_rmats_silrest2['start3ss'] + window

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & (2*window > (df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss'])) & \
                    ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) > window),
                    'start'] = df_rmats_silrest2['start3ss'] - window
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & (2*window > (df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss'])) & \
                    ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) > window),
                    'end'] = df_rmats_silrest2['start5ss'] + window

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) < window),
                    'start'] = df_rmats_silrest2['start3ss']
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) < window),
                    'end'] = df_rmats_silrest2['start5ss']

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) < window),
                    'start'] = df_rmats_silrest2['start5ss']
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) < window),
                    'end'] = df_rmats_silrest2['start3ss']

    df_rmats_silrest2['len'] = df_rmats_silrest2['end'] - df_rmats_silrest2['start']

    df_rmats_silrest2['diff'] = abs(df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss'])

    df_rmats_silrest2 = df_rmats_silrest2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_silrest2 = df_rmats_silrest2[~(df_rmats_silrest2['start'].isna() | df_rmats_silrest2['end'].isna())]

    df_rmats_silrest2['start'] = df_rmats_silrest2['start'].astype(int)
    df_rmats_silrest2['end'] = df_rmats_silrest2['end'].astype(int)

    rmats_silrest2 = pbt.BedTool.from_dataframe(df_rmats_silrest2)
    rmats_silrest2 = rmats_silrest2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_silrest2 = rmats_silrest2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_silrest2['base_coverage_5ss']
    del df_rmats_silrest2['chrom2']
    del df_rmats_silrest2['start2']
    del df_rmats_silrest2['end2']
    del df_rmats_silrest2['useless']
    del df_rmats_silrest2['strand2']

    df_rmats_silrest2_grouped10 = df_rmats_silrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_silrest2_grouped12 = df_rmats_silrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_silrest2 = df_rmats_silrest2_grouped10.join(df_rmats_silrest2_grouped12).reset_index()
    df_rmats_silrest2 = df_rmats_silrest2.set_index('index')

    df_rmats_silrest_overlap = df_rmats_silrest.join(df_rmats_silrest2)

    df_rmats_silrest_overlap = df_rmats_silrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_silrest_overlap.loc[df_rmats_silrest_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_silrest_overlap.loc[df_rmats_silrest_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_silrest_overlap.loc[df_rmats_silrest_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_silrest_overlap['cDNA_coverage_5ss'] = df_rmats_silrest_overlap['cDNA_coverage_5ss'] - df_rmats_silrest_overlap['overlap_score']
    df_rmats_silrest_overlap['cDNA_coverage_3ss'] = df_rmats_silrest_overlap['cDNA_coverage_3ss'] - df_rmats_silrest_overlap['overlap_score']
    df_rmats_silrest_overlap['base_coverage_5ss'] = df_rmats_silrest_overlap['base_coverage_5ss'] - df_rmats_silrest_overlap['overlap_base']
    df_rmats_silrest_overlap['base_coverage_3ss'] = df_rmats_silrest_overlap['base_coverage_3ss'] - df_rmats_silrest_overlap['overlap_base']

    df_rmats_silrest = df_rmats_silrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss.reset_index()
    rmats_ctrl_3ss = pbt.BedTool.from_dataframe(df_rmats_ctrl_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_ctrl_3ss = rmats_ctrl_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_ctrl_3ss = rmats_ctrl_3ss.to_dataframe(header=None)
    df_ctrl_3ss = df_ctrl_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_ctrl_3ss_grouped10 = df_ctrl_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_ctrl_3ss_grouped_12 = df_ctrl_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_ctrl_3ss = df_ctrl_3ss_grouped10.join(df_ctrl_3ss_grouped_12).reset_index()
    df_ctrl_3ss = df_ctrl_3ss.set_index(3)
    df_rmats_ctrl_3ss_temp = pd.merge(df_rmats_ctrl_3ss_temp, df_ctrl_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_ctrl_3ss_temp['index']
    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss_temp.set_index(14)
    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss.reset_index()
    rmats_ctrl_5ss = pbt.BedTool.from_dataframe(df_rmats_ctrl_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_ctrl_5ss = rmats_ctrl_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_ctrl_5ss = rmats_ctrl_5ss.to_dataframe(header=None)
    df_ctrl_5ss = df_ctrl_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_ctrl_5ss_grouped10 = df_ctrl_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_ctrl_5ss_grouped_12 = df_ctrl_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_ctrl_5ss = df_ctrl_5ss_grouped10.join(df_ctrl_5ss_grouped_12).reset_index()
    df_ctrl_5ss = df_ctrl_5ss.set_index(3)
    df_rmats_ctrl_5ss_temp = pd.merge(df_rmats_ctrl_5ss_temp, df_ctrl_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_ctrl_5ss_temp['index']
    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss_temp.set_index(14)
    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_ctrl = df_rmats_ctrl_5ss_temp.join(df_rmats_ctrl_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_ctrl2 = df_rmats_ctrl[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_ctrl2 = df_rmats_ctrl2.reset_index()

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & (2*window > (df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss'])) & \
                    ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) > window),
                    'start'] = df_rmats_ctrl2['start5ss'] - window
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & (2*window > (df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss'])) & \
                    ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) > window),
                    'end'] = df_rmats_ctrl2['start3ss'] + window

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & (2*window > (df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss'])) & \
                    ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) > window),
                    'start'] = df_rmats_ctrl2['start3ss'] - window
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & (2*window > (df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss'])) & \
                    ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) > window),
                    'end'] = df_rmats_ctrl2['start5ss'] + window

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) < window),
                    'start'] = df_rmats_ctrl2['start3ss']
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) < window),
                    'end'] = df_rmats_ctrl2['start5ss']

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) < window),
                    'start'] = df_rmats_ctrl2['start5ss']
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) < window),
                    'end'] = df_rmats_ctrl2['start3ss']

    df_rmats_ctrl2['len'] = df_rmats_ctrl2['end'] - df_rmats_ctrl2['start']

    df_rmats_ctrl2['diff'] = abs(df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss'])

    df_rmats_ctrl2 = df_rmats_ctrl2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_ctrl2 = df_rmats_ctrl2[~(df_rmats_ctrl2['start'].isna() | df_rmats_ctrl2['end'].isna())]

    df_rmats_ctrl2['start'] = df_rmats_ctrl2['start'].astype(int)
    df_rmats_ctrl2['end'] = df_rmats_ctrl2['end'].astype(int)

    rmats_ctrl2 = pbt.BedTool.from_dataframe(df_rmats_ctrl2)
    rmats_ctrl2 = rmats_ctrl2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_ctrl2 = rmats_ctrl2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_ctrl2['base_coverage_5ss']
    del df_rmats_ctrl2['chrom2']
    del df_rmats_ctrl2['start2']
    del df_rmats_ctrl2['end2']
    del df_rmats_ctrl2['useless']
    del df_rmats_ctrl2['strand2']

    df_rmats_ctrl2_grouped10 = df_rmats_ctrl2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_ctrl2_grouped12 = df_rmats_ctrl2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_ctrl2 = df_rmats_ctrl2_grouped10.join(df_rmats_ctrl2_grouped12).reset_index()
    df_rmats_ctrl2 = df_rmats_ctrl2.set_index('index')

    df_rmats_ctrl_overlap = df_rmats_ctrl.join(df_rmats_ctrl2)

    df_rmats_ctrl_overlap = df_rmats_ctrl_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_ctrl_overlap.loc[df_rmats_ctrl_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_ctrl_overlap.loc[df_rmats_ctrl_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_ctrl_overlap.loc[df_rmats_ctrl_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_ctrl_overlap['cDNA_coverage_5ss'] = df_rmats_ctrl_overlap['cDNA_coverage_5ss'] - df_rmats_ctrl_overlap['overlap_score']
    df_rmats_ctrl_overlap['cDNA_coverage_3ss'] = df_rmats_ctrl_overlap['cDNA_coverage_3ss'] - df_rmats_ctrl_overlap['overlap_score']
    df_rmats_ctrl_overlap['base_coverage_5ss'] = df_rmats_ctrl_overlap['base_coverage_5ss'] - df_rmats_ctrl_overlap['overlap_base']
    df_rmats_ctrl_overlap['base_coverage_3ss'] = df_rmats_ctrl_overlap['base_coverage_3ss'] - df_rmats_ctrl_overlap['overlap_base']

    df_rmats_ctrl = df_rmats_ctrl_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_const_3ss_temp = df_rmats_const_3ss.reset_index()
    rmats_const_3ss = pbt.BedTool.from_dataframe(df_rmats_const_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_const_3ss = rmats_const_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_const_3ss = rmats_const_3ss.to_dataframe(header=None)
    df_const_3ss = df_const_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_const_3ss_grouped10 = df_const_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_const_3ss_grouped_12 = df_const_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_const_3ss = df_const_3ss_grouped10.join(df_const_3ss_grouped_12).reset_index()
    df_const_3ss = df_const_3ss.set_index(3)
    df_rmats_const_3ss_temp = pd.merge(df_rmats_const_3ss_temp, df_const_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_const_3ss_temp['index']
    df_rmats_const_3ss_temp = df_rmats_const_3ss_temp.set_index(14)
    df_rmats_const_3ss_temp = df_rmats_const_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_const_5ss_temp = df_rmats_const_5ss.reset_index()
    rmats_const_5ss = pbt.BedTool.from_dataframe(df_rmats_const_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_const_5ss = rmats_const_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_const_5ss = rmats_const_5ss.to_dataframe(header=None)
    df_const_5ss = df_const_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_const_5ss_grouped10 = df_const_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_const_5ss_grouped_12 = df_const_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_const_5ss = df_const_5ss_grouped10.join(df_const_5ss_grouped_12).reset_index()
    df_const_5ss = df_const_5ss.set_index(3)
    df_rmats_const_5ss_temp = pd.merge(df_rmats_const_5ss_temp, df_const_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_const_5ss_temp['index']
    df_rmats_const_5ss_temp = df_rmats_const_5ss_temp.set_index(14)
    df_rmats_const_5ss_temp = df_rmats_const_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_const = df_rmats_const_5ss_temp.join(df_rmats_const_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_const2 = df_rmats_const[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_const2 = df_rmats_const2.reset_index()

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & (2*window > (df_rmats_const2['start5ss'] - df_rmats_const2['start3ss'])) & \
                    ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) > window),
                    'start'] = df_rmats_const2['start5ss'] - window
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & (2*window > (df_rmats_const2['start5ss'] - df_rmats_const2['start3ss'])) & \
                    ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) > window),
                    'end'] = df_rmats_const2['start3ss'] + window

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & (2*window > (df_rmats_const2['start3ss'] - df_rmats_const2['start5ss'])) & \
                    ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) > window),
                    'start'] = df_rmats_const2['start3ss'] - window
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & (2*window > (df_rmats_const2['start3ss'] - df_rmats_const2['start5ss'])) & \
                    ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) > window),
                    'end'] = df_rmats_const2['start5ss'] + window

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) < window),
                    'start'] = df_rmats_const2['start3ss']
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) < window),
                    'end'] = df_rmats_const2['start5ss']

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) < window),
                    'start'] = df_rmats_const2['start5ss']
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) < window),
                    'end'] = df_rmats_const2['start3ss']

    df_rmats_const2['len'] = df_rmats_const2['end'] - df_rmats_const2['start']

    df_rmats_const2['diff'] = abs(df_rmats_const2['start3ss'] - df_rmats_const2['start5ss'])

    df_rmats_const2 = df_rmats_const2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_const2 = df_rmats_const2[~(df_rmats_const2['start'].isna() | df_rmats_const2['end'].isna())]

    df_rmats_const2['start'] = df_rmats_const2['start'].astype(int)
    df_rmats_const2['end'] = df_rmats_const2['end'].astype(int)

    rmats_const2 = pbt.BedTool.from_dataframe(df_rmats_const2)
    rmats_const2 = rmats_const2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_const2 = rmats_const2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_const2['base_coverage_5ss']
    del df_rmats_const2['chrom2']
    del df_rmats_const2['start2']
    del df_rmats_const2['end2']
    del df_rmats_const2['useless']
    del df_rmats_const2['strand2']

    df_rmats_const2_grouped10 = df_rmats_const2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_const2_grouped12 = df_rmats_const2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_const2 = df_rmats_const2_grouped10.join(df_rmats_const2_grouped12).reset_index()
    df_rmats_const2 = df_rmats_const2.set_index('index')

    df_rmats_const_overlap = df_rmats_const.join(df_rmats_const2)

    df_rmats_const_overlap = df_rmats_const_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_const_overlap.loc[df_rmats_const_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_const_overlap.loc[df_rmats_const_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_const_overlap.loc[df_rmats_const_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_const_overlap['cDNA_coverage_5ss'] = df_rmats_const_overlap['cDNA_coverage_5ss'] - df_rmats_const_overlap['overlap_score']
    df_rmats_const_overlap['cDNA_coverage_3ss'] = df_rmats_const_overlap['cDNA_coverage_3ss'] - df_rmats_const_overlap['overlap_score']
    df_rmats_const_overlap['base_coverage_5ss'] = df_rmats_const_overlap['base_coverage_5ss'] - df_rmats_const_overlap['overlap_base']
    df_rmats_const_overlap['base_coverage_3ss'] = df_rmats_const_overlap['base_coverage_3ss'] - df_rmats_const_overlap['overlap_base']

    df_rmats_const = df_rmats_const_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    
    df_rmats_enh_3ss_temp = df_rmats_enh_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss.reset_index()
    rmats_enh_5ss = pbt.BedTool.from_dataframe(df_rmats_enh_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enh_5ss = rmats_enh_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enh_5ss = rmats_enh_5ss.to_dataframe(header=None)
    df_enh_5ss = df_enh_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enh_5ss_grouped10 = df_enh_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enh_5ss_grouped_12 = df_enh_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enh_5ss = df_enh_5ss_grouped10.join(df_enh_5ss_grouped_12).reset_index()
    df_enh_5ss = df_enh_5ss.set_index(3)
    df_rmats_enh_5ss_temp = pd.merge(df_rmats_enh_5ss_temp, df_enh_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enh_5ss_temp['index']
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss_temp.set_index(14)
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_enh = df_rmats_enh_5ss_temp.join(df_rmats_enh_3ss_temp, lsuffix='5ss', rsuffix='3ss')


    df_rmats_enh2 = df_rmats_enh[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
       'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_enh2 = df_rmats_enh2.reset_index()

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & (2*window > (df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss'])) & \
                      ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) > window),
                      'start'] = df_rmats_enh2['start5ss'] - window
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & (2*window > (df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss'])) & \
                      ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) > window),
                      'end'] = df_rmats_enh2['start3ss'] + window

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & (2*window > (df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss'])) & \
                      ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) > window),
                      'start'] = df_rmats_enh2['start3ss'] - window
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & (2*window > (df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss'])) & \
                      ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) > window),
                      'end'] = df_rmats_enh2['start5ss'] + window

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) < window),
                      'start'] = df_rmats_enh2['start3ss']
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '+') & ((df_rmats_enh2['start5ss'] - df_rmats_enh2['start3ss']) < window),
                      'end'] = df_rmats_enh2['start5ss']

    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) < window),
                      'start'] = df_rmats_enh2['start5ss']
    df_rmats_enh2.loc[(df_rmats_enh2['strand5ss'] == '-') & ((df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss']) < window),
                      'end'] = df_rmats_enh2['start3ss']

    df_rmats_enh2['len'] = df_rmats_enh2['end'] - df_rmats_enh2['start']

    df_rmats_enh2['diff'] = abs(df_rmats_enh2['start3ss'] - df_rmats_enh2['start5ss'])

    df_rmats_enh2 = df_rmats_enh2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_enh2 = df_rmats_enh2[~(df_rmats_enh2['start'].isna() | df_rmats_enh2['end'].isna())]

    df_rmats_enh2['start'] = df_rmats_enh2['start'].astype(int)
    df_rmats_enh2['end'] = df_rmats_enh2['end'].astype(int)

    rmats_enh2 = pbt.BedTool.from_dataframe(df_rmats_enh2)
    rmats_enh2 = rmats_enh2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_enh2 = rmats_enh2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_enh2['base_coverage_5ss']
    del df_rmats_enh2['chrom2']
    del df_rmats_enh2['start2']
    del df_rmats_enh2['end2']
    del df_rmats_enh2['useless']
    del df_rmats_enh2['strand2']

    df_rmats_enh2_grouped10 = df_rmats_enh2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_enh2_grouped12 = df_rmats_enh2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_enh2 = df_rmats_enh2_grouped10.join(df_rmats_enh2_grouped12).reset_index()
    df_rmats_enh2 = df_rmats_enh2.set_index('index')

    df_rmats_enh_overlap = df_rmats_enh.join(df_rmats_enh2)

    df_rmats_enh_overlap = df_rmats_enh_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_enh_overlap.loc[df_rmats_enh_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_enh_overlap.loc[df_rmats_enh_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_enh_overlap.loc[df_rmats_enh_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_enh_overlap['cDNA_coverage_5ss'] = df_rmats_enh_overlap['cDNA_coverage_5ss'] - df_rmats_enh_overlap['overlap_score']
    df_rmats_enh_overlap['cDNA_coverage_3ss'] = df_rmats_enh_overlap['cDNA_coverage_3ss'] - df_rmats_enh_overlap['overlap_score']
    df_rmats_enh_overlap['base_coverage_5ss'] = df_rmats_enh_overlap['base_coverage_5ss'] - df_rmats_enh_overlap['overlap_base']
    df_rmats_enh_overlap['base_coverage_3ss'] = df_rmats_enh_overlap['base_coverage_3ss'] - df_rmats_enh_overlap['overlap_base']

    df_rmats_enh = df_rmats_enh_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_sil_3ss_temp = df_rmats_sil_3ss.reset_index()
    rmats_sil_3ss = pbt.BedTool.from_dataframe(df_rmats_sil_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_sil_3ss = rmats_sil_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_sil_3ss = rmats_sil_3ss.to_dataframe(header=None)
    df_sil_3ss = df_sil_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_sil_3ss_grouped10 = df_sil_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_sil_3ss_grouped_12 = df_sil_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_sil_3ss = df_sil_3ss_grouped10.join(df_sil_3ss_grouped_12).reset_index()
    df_sil_3ss = df_sil_3ss.set_index(3)
    df_rmats_sil_3ss_temp = pd.merge(df_rmats_sil_3ss_temp, df_sil_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_sil_3ss_temp['index']
    df_rmats_sil_3ss_temp = df_rmats_sil_3ss_temp.set_index(14)
    df_rmats_sil_3ss_temp = df_rmats_sil_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_sil_5ss_temp = df_rmats_sil_5ss.reset_index()
    rmats_sil_5ss = pbt.BedTool.from_dataframe(df_rmats_sil_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_sil_5ss = rmats_sil_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_sil_5ss = rmats_sil_5ss.to_dataframe(header=None)
    df_sil_5ss = df_sil_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_sil_5ss_grouped10 = df_sil_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_sil_5ss_grouped_12 = df_sil_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_sil_5ss = df_sil_5ss_grouped10.join(df_sil_5ss_grouped_12).reset_index()
    df_sil_5ss = df_sil_5ss.set_index(3)
    df_rmats_sil_5ss_temp = pd.merge(df_rmats_sil_5ss_temp, df_sil_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_sil_5ss_temp['index']
    df_rmats_sil_5ss_temp = df_rmats_sil_5ss_temp.set_index(14)
    df_rmats_sil_5ss_temp = df_rmats_sil_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_sil = df_rmats_sil_5ss_temp.join(df_rmats_sil_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_sil2 = df_rmats_sil[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
       'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_sil2 = df_rmats_sil2.reset_index()

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & (2*window > (df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss'])) & \
                      ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) > window),
                      'start'] = df_rmats_sil2['start5ss'] - window
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & (2*window > (df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss'])) & \
                      ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) > window),
                      'end'] = df_rmats_sil2['start3ss'] + window

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & (2*window > (df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss'])) & \
                      ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) > window),
                      'start'] = df_rmats_sil2['start3ss'] - window
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & (2*window > (df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss'])) & \
                      ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) > window),
                      'end'] = df_rmats_sil2['start5ss'] + window

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) < window),
                      'start'] = df_rmats_sil2['start3ss']
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '+') & ((df_rmats_sil2['start5ss'] - df_rmats_sil2['start3ss']) < window),
                      'end'] = df_rmats_sil2['start5ss']

    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) < window),
                      'start'] = df_rmats_sil2['start5ss']
    df_rmats_sil2.loc[(df_rmats_sil2['strand5ss'] == '-') & ((df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss']) < window),
                      'end'] = df_rmats_sil2['start3ss']

    df_rmats_sil2['len'] = df_rmats_sil2['end'] - df_rmats_sil2['start']

    df_rmats_sil2['diff'] = abs(df_rmats_sil2['start3ss'] - df_rmats_sil2['start5ss'])

    df_rmats_sil2 = df_rmats_sil2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_sil2 = df_rmats_sil2[~(df_rmats_sil2['start'].isna() | df_rmats_sil2['end'].isna())]

    df_rmats_sil2['start'] = df_rmats_sil2['start'].astype(int)
    df_rmats_sil2['end'] = df_rmats_sil2['end'].astype(int)

    rmats_sil2 = pbt.BedTool.from_dataframe(df_rmats_sil2)
    rmats_sil2 = rmats_sil2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_sil2 = rmats_sil2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_sil2['base_coverage_5ss']
    del df_rmats_sil2['chrom2']
    del df_rmats_sil2['start2']
    del df_rmats_sil2['end2']
    del df_rmats_sil2['useless']
    del df_rmats_sil2['strand2']

    df_rmats_sil2_grouped10 = df_rmats_sil2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_sil2_grouped12 = df_rmats_sil2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_sil2 = df_rmats_sil2_grouped10.join(df_rmats_sil2_grouped12).reset_index()
    df_rmats_sil2 = df_rmats_sil2.set_index('index')

    df_rmats_sil_overlap = df_rmats_sil.join(df_rmats_sil2)

    df_rmats_sil_overlap = df_rmats_sil_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_sil_overlap.loc[df_rmats_sil_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_sil_overlap.loc[df_rmats_sil_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_sil_overlap.loc[df_rmats_sil_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_sil_overlap['cDNA_coverage_5ss'] = df_rmats_sil_overlap['cDNA_coverage_5ss'] - df_rmats_sil_overlap['overlap_score']
    df_rmats_sil_overlap['cDNA_coverage_3ss'] = df_rmats_sil_overlap['cDNA_coverage_3ss'] - df_rmats_sil_overlap['overlap_score']
    df_rmats_sil_overlap['base_coverage_5ss'] = df_rmats_sil_overlap['base_coverage_5ss'] - df_rmats_sil_overlap['overlap_base']
    df_rmats_sil_overlap['base_coverage_3ss'] = df_rmats_sil_overlap['base_coverage_3ss'] - df_rmats_sil_overlap['overlap_base']

    df_rmats_sil = df_rmats_sil_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss.reset_index()
    rmats_enhrest_3ss = pbt.BedTool.from_dataframe(df_rmats_enhrest_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enhrest_3ss = rmats_enhrest_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enhrest_3ss = rmats_enhrest_3ss.to_dataframe(header=None)
    df_enhrest_3ss = df_enhrest_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enhrest_3ss_grouped10 = df_enhrest_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enhrest_3ss_grouped_12 = df_enhrest_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enhrest_3ss = df_enhrest_3ss_grouped10.join(df_enhrest_3ss_grouped_12).reset_index()
    df_enhrest_3ss = df_enhrest_3ss.set_index(3)
    df_rmats_enhrest_3ss_temp = pd.merge(df_rmats_enhrest_3ss_temp, df_enhrest_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enhrest_3ss_temp['index']
    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss_temp.set_index(14)
    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss.reset_index()
    rmats_enhrest_5ss = pbt.BedTool.from_dataframe(df_rmats_enhrest_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_enhrest_5ss = rmats_enhrest_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_enhrest_5ss = rmats_enhrest_5ss.to_dataframe(header=None)
    df_enhrest_5ss = df_enhrest_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_enhrest_5ss_grouped10 = df_enhrest_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_enhrest_5ss_grouped_12 = df_enhrest_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_enhrest_5ss = df_enhrest_5ss_grouped10.join(df_enhrest_5ss_grouped_12).reset_index()
    df_enhrest_5ss = df_enhrest_5ss.set_index(3)
    df_rmats_enhrest_5ss_temp = pd.merge(df_rmats_enhrest_5ss_temp, df_enhrest_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_enhrest_5ss_temp['index']
    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss_temp.set_index(14)
    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_enhrest = df_rmats_enhrest_5ss_temp.join(df_rmats_enhrest_3ss_temp, lsuffix='5ss', rsuffix='3ss')

    df_rmats_enhrest2 = df_rmats_enhrest[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_enhrest2 = df_rmats_enhrest2.reset_index()

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & (2*window > (df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss'])) & \
                    ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) > window),
                    'start'] = df_rmats_enhrest2['start5ss'] - window
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & (2*window > (df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss'])) & \
                    ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) > window),
                    'end'] = df_rmats_enhrest2['start3ss'] + window

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & (2*window > (df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss'])) & \
                    ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) > window),
                    'start'] = df_rmats_enhrest2['start3ss'] - window
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & (2*window > (df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss'])) & \
                    ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) > window),
                    'end'] = df_rmats_enhrest2['start5ss'] + window

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) < window),
                    'start'] = df_rmats_enhrest2['start3ss']
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '+') & ((df_rmats_enhrest2['start5ss'] - df_rmats_enhrest2['start3ss']) < window),
                    'end'] = df_rmats_enhrest2['start5ss']

    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) < window),
                    'start'] = df_rmats_enhrest2['start5ss']
    df_rmats_enhrest2.loc[(df_rmats_enhrest2['strand5ss'] == '-') & ((df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss']) < window),
                    'end'] = df_rmats_enhrest2['start3ss']

    df_rmats_enhrest2['len'] = df_rmats_enhrest2['end'] - df_rmats_enhrest2['start']

    df_rmats_enhrest2['diff'] = abs(df_rmats_enhrest2['start3ss'] - df_rmats_enhrest2['start5ss'])

    df_rmats_enhrest2 = df_rmats_enhrest2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_enhrest2 = df_rmats_enhrest2[~(df_rmats_enhrest2['start'].isna() | df_rmats_enhrest2['end'].isna())]

    df_rmats_enhrest2['start'] = df_rmats_enhrest2['start'].astype(int)
    df_rmats_enhrest2['end'] = df_rmats_enhrest2['end'].astype(int)

    rmats_enhrest2 = pbt.BedTool.from_dataframe(df_rmats_enhrest2)
    rmats_enhrest2 = rmats_enhrest2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_enhrest2 = rmats_enhrest2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_enhrest2['base_coverage_5ss']
    del df_rmats_enhrest2['chrom2']
    del df_rmats_enhrest2['start2']
    del df_rmats_enhrest2['end2']
    del df_rmats_enhrest2['useless']
    del df_rmats_enhrest2['strand2']

    df_rmats_enhrest2_grouped10 = df_rmats_enhrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_enhrest2_grouped12 = df_rmats_enhrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_enhrest2 = df_rmats_enhrest2_grouped10.join(df_rmats_enhrest2_grouped12).reset_index()
    df_rmats_enhrest2 = df_rmats_enhrest2.set_index('index')

    df_rmats_enhrest_overlap = df_rmats_enhrest.join(df_rmats_enhrest2)

    df_rmats_enhrest_overlap = df_rmats_enhrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_enhrest_overlap.loc[df_rmats_enhrest_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_enhrest_overlap.loc[df_rmats_enhrest_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_enhrest_overlap.loc[df_rmats_enhrest_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_enhrest_overlap['cDNA_coverage_5ss'] = df_rmats_enhrest_overlap['cDNA_coverage_5ss'] - df_rmats_enhrest_overlap['overlap_score']
    df_rmats_enhrest_overlap['cDNA_coverage_3ss'] = df_rmats_enhrest_overlap['cDNA_coverage_3ss'] - df_rmats_enhrest_overlap['overlap_score']
    df_rmats_enhrest_overlap['base_coverage_5ss'] = df_rmats_enhrest_overlap['base_coverage_5ss'] - df_rmats_enhrest_overlap['overlap_base']
    df_rmats_enhrest_overlap['base_coverage_3ss'] = df_rmats_enhrest_overlap['base_coverage_3ss'] - df_rmats_enhrest_overlap['overlap_base']

    df_rmats_enhrest = df_rmats_enhrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    
    
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss.reset_index()
    rmats_silrest_3ss = pbt.BedTool.from_dataframe(df_rmats_silrest_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_silrest_3ss = rmats_silrest_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_silrest_3ss = rmats_silrest_3ss.to_dataframe(header=None)
    df_silrest_3ss = df_silrest_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_silrest_3ss_grouped10 = df_silrest_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_silrest_3ss_grouped_12 = df_silrest_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_silrest_3ss = df_silrest_3ss_grouped10.join(df_silrest_3ss_grouped_12).reset_index()
    df_silrest_3ss = df_silrest_3ss.set_index(3)
    df_rmats_silrest_3ss_temp = pd.merge(df_rmats_silrest_3ss_temp, df_silrest_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_silrest_3ss_temp['index']
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss_temp.set_index(14)
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss.reset_index()
    rmats_silrest_5ss = pbt.BedTool.from_dataframe(df_rmats_silrest_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_silrest_5ss = rmats_silrest_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_silrest_5ss = rmats_silrest_5ss.to_dataframe(header=None)
    df_silrest_5ss = df_silrest_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_silrest_5ss_grouped10 = df_silrest_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_silrest_5ss_grouped_12 = df_silrest_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_silrest_5ss = df_silrest_5ss_grouped10.join(df_silrest_5ss_grouped_12).reset_index()
    df_silrest_5ss = df_silrest_5ss.set_index(3)
    df_rmats_silrest_5ss_temp = pd.merge(df_rmats_silrest_5ss_temp, df_silrest_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_silrest_5ss_temp['index']
    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss_temp.set_index(14)
    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_silrest = df_rmats_silrest_5ss_temp.join(df_rmats_silrest_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_silrest2 = df_rmats_silrest[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_silrest2 = df_rmats_silrest2.reset_index()

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & (2*window > (df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss'])) & \
                    ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) > window),
                    'start'] = df_rmats_silrest2['start5ss'] - window
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & (2*window > (df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss'])) & \
                    ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) > window),
                    'end'] = df_rmats_silrest2['start3ss'] + window

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & (2*window > (df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss'])) & \
                    ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) > window),
                    'start'] = df_rmats_silrest2['start3ss'] - window
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & (2*window > (df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss'])) & \
                    ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) > window),
                    'end'] = df_rmats_silrest2['start5ss'] + window

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) < window),
                    'start'] = df_rmats_silrest2['start3ss']
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '+') & ((df_rmats_silrest2['start5ss'] - df_rmats_silrest2['start3ss']) < window),
                    'end'] = df_rmats_silrest2['start5ss']

    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) < window),
                    'start'] = df_rmats_silrest2['start5ss']
    df_rmats_silrest2.loc[(df_rmats_silrest2['strand5ss'] == '-') & ((df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss']) < window),
                    'end'] = df_rmats_silrest2['start3ss']

    df_rmats_silrest2['len'] = df_rmats_silrest2['end'] - df_rmats_silrest2['start']

    df_rmats_silrest2['diff'] = abs(df_rmats_silrest2['start3ss'] - df_rmats_silrest2['start5ss'])

    df_rmats_silrest2 = df_rmats_silrest2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_silrest2 = df_rmats_silrest2[~(df_rmats_silrest2['start'].isna() | df_rmats_silrest2['end'].isna())]

    df_rmats_silrest2['start'] = df_rmats_silrest2['start'].astype(int)
    df_rmats_silrest2['end'] = df_rmats_silrest2['end'].astype(int)

    rmats_silrest2 = pbt.BedTool.from_dataframe(df_rmats_silrest2)
    rmats_silrest2 = rmats_silrest2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_silrest2 = rmats_silrest2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_silrest2['base_coverage_5ss']
    del df_rmats_silrest2['chrom2']
    del df_rmats_silrest2['start2']
    del df_rmats_silrest2['end2']
    del df_rmats_silrest2['useless']
    del df_rmats_silrest2['strand2']

    df_rmats_silrest2_grouped10 = df_rmats_silrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_silrest2_grouped12 = df_rmats_silrest2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_silrest2 = df_rmats_silrest2_grouped10.join(df_rmats_silrest2_grouped12).reset_index()
    df_rmats_silrest2 = df_rmats_silrest2.set_index('index')

    df_rmats_silrest_overlap = df_rmats_silrest.join(df_rmats_silrest2)

    df_rmats_silrest_overlap = df_rmats_silrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_silrest_overlap.loc[df_rmats_silrest_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_silrest_overlap.loc[df_rmats_silrest_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_silrest_overlap.loc[df_rmats_silrest_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_silrest_overlap['cDNA_coverage_5ss'] = df_rmats_silrest_overlap['cDNA_coverage_5ss'] - df_rmats_silrest_overlap['overlap_score']
    df_rmats_silrest_overlap['cDNA_coverage_3ss'] = df_rmats_silrest_overlap['cDNA_coverage_3ss'] - df_rmats_silrest_overlap['overlap_score']
    df_rmats_silrest_overlap['base_coverage_5ss'] = df_rmats_silrest_overlap['base_coverage_5ss'] - df_rmats_silrest_overlap['overlap_base']
    df_rmats_silrest_overlap['base_coverage_3ss'] = df_rmats_silrest_overlap['base_coverage_3ss'] - df_rmats_silrest_overlap['overlap_base']

    df_rmats_silrest = df_rmats_silrest_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss.reset_index()
    rmats_ctrl_3ss = pbt.BedTool.from_dataframe(df_rmats_ctrl_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_ctrl_3ss = rmats_ctrl_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_ctrl_3ss = rmats_ctrl_3ss.to_dataframe(header=None)
    df_ctrl_3ss = df_ctrl_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_ctrl_3ss_grouped10 = df_ctrl_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_ctrl_3ss_grouped_12 = df_ctrl_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_ctrl_3ss = df_ctrl_3ss_grouped10.join(df_ctrl_3ss_grouped_12).reset_index()
    df_ctrl_3ss = df_ctrl_3ss.set_index(3)
    df_rmats_ctrl_3ss_temp = pd.merge(df_rmats_ctrl_3ss_temp, df_ctrl_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_ctrl_3ss_temp['index']
    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss_temp.set_index(14)
    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss.reset_index()
    rmats_ctrl_5ss = pbt.BedTool.from_dataframe(df_rmats_ctrl_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_ctrl_5ss = rmats_ctrl_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_ctrl_5ss = rmats_ctrl_5ss.to_dataframe(header=None)
    df_ctrl_5ss = df_ctrl_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_ctrl_5ss_grouped10 = df_ctrl_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_ctrl_5ss_grouped_12 = df_ctrl_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_ctrl_5ss = df_ctrl_5ss_grouped10.join(df_ctrl_5ss_grouped_12).reset_index()
    df_ctrl_5ss = df_ctrl_5ss.set_index(3)
    df_rmats_ctrl_5ss_temp = pd.merge(df_rmats_ctrl_5ss_temp, df_ctrl_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_ctrl_5ss_temp['index']
    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss_temp.set_index(14)
    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_ctrl = df_rmats_ctrl_5ss_temp.join(df_rmats_ctrl_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_ctrl2 = df_rmats_ctrl[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_ctrl2 = df_rmats_ctrl2.reset_index()

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & (2*window > (df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss'])) & \
                    ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) > window),
                    'start'] = df_rmats_ctrl2['start5ss'] - window
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & (2*window > (df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss'])) & \
                    ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) > window),
                    'end'] = df_rmats_ctrl2['start3ss'] + window

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & (2*window > (df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss'])) & \
                    ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) > window),
                    'start'] = df_rmats_ctrl2['start3ss'] - window
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & (2*window > (df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss'])) & \
                    ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) > window),
                    'end'] = df_rmats_ctrl2['start5ss'] + window

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) < window),
                    'start'] = df_rmats_ctrl2['start3ss']
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '+') & ((df_rmats_ctrl2['start5ss'] - df_rmats_ctrl2['start3ss']) < window),
                    'end'] = df_rmats_ctrl2['start5ss']

    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) < window),
                    'start'] = df_rmats_ctrl2['start5ss']
    df_rmats_ctrl2.loc[(df_rmats_ctrl2['strand5ss'] == '-') & ((df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss']) < window),
                    'end'] = df_rmats_ctrl2['start3ss']

    df_rmats_ctrl2['len'] = df_rmats_ctrl2['end'] - df_rmats_ctrl2['start']

    df_rmats_ctrl2['diff'] = abs(df_rmats_ctrl2['start3ss'] - df_rmats_ctrl2['start5ss'])

    df_rmats_ctrl2 = df_rmats_ctrl2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_ctrl2 = df_rmats_ctrl2[~(df_rmats_ctrl2['start'].isna() | df_rmats_ctrl2['end'].isna())]

    df_rmats_ctrl2['start'] = df_rmats_ctrl2['start'].astype(int)
    df_rmats_ctrl2['end'] = df_rmats_ctrl2['end'].astype(int)

    rmats_ctrl2 = pbt.BedTool.from_dataframe(df_rmats_ctrl2)
    rmats_ctrl2 = rmats_ctrl2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_ctrl2 = rmats_ctrl2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_ctrl2['base_coverage_5ss']
    del df_rmats_ctrl2['chrom2']
    del df_rmats_ctrl2['start2']
    del df_rmats_ctrl2['end2']
    del df_rmats_ctrl2['useless']
    del df_rmats_ctrl2['strand2']

    df_rmats_ctrl2_grouped10 = df_rmats_ctrl2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_ctrl2_grouped12 = df_rmats_ctrl2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_ctrl2 = df_rmats_ctrl2_grouped10.join(df_rmats_ctrl2_grouped12).reset_index()
    df_rmats_ctrl2 = df_rmats_ctrl2.set_index('index')

    df_rmats_ctrl_overlap = df_rmats_ctrl.join(df_rmats_ctrl2)

    df_rmats_ctrl_overlap = df_rmats_ctrl_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_ctrl_overlap.loc[df_rmats_ctrl_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_ctrl_overlap.loc[df_rmats_ctrl_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_ctrl_overlap.loc[df_rmats_ctrl_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_ctrl_overlap['cDNA_coverage_5ss'] = df_rmats_ctrl_overlap['cDNA_coverage_5ss'] - df_rmats_ctrl_overlap['overlap_score']
    df_rmats_ctrl_overlap['cDNA_coverage_3ss'] = df_rmats_ctrl_overlap['cDNA_coverage_3ss'] - df_rmats_ctrl_overlap['overlap_score']
    df_rmats_ctrl_overlap['base_coverage_5ss'] = df_rmats_ctrl_overlap['base_coverage_5ss'] - df_rmats_ctrl_overlap['overlap_base']
    df_rmats_ctrl_overlap['base_coverage_3ss'] = df_rmats_ctrl_overlap['base_coverage_3ss'] - df_rmats_ctrl_overlap['overlap_base']

    df_rmats_ctrl = df_rmats_ctrl_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_const_3ss_temp = df_rmats_const_3ss.reset_index()
    rmats_const_3ss = pbt.BedTool.from_dataframe(df_rmats_const_3ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_const_3ss = rmats_const_3ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_const_3ss = rmats_const_3ss.to_dataframe(header=None)
    df_const_3ss = df_const_3ss[[0, 1, 2, 3, 4, 10, 12]]
    df_const_3ss_grouped10 = df_const_3ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_const_3ss_grouped_12 = df_const_3ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_const_3ss = df_const_3ss_grouped10.join(df_const_3ss_grouped_12).reset_index()
    df_const_3ss = df_const_3ss.set_index(3)
    df_rmats_const_3ss_temp = pd.merge(df_rmats_const_3ss_temp, df_const_3ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_const_3ss_temp['index']
    df_rmats_const_3ss_temp = df_rmats_const_3ss_temp.set_index(14)
    df_rmats_const_3ss_temp = df_rmats_const_3ss_temp.rename(columns={10: 'cDNA_coverage_3ss', 12: 'base_coverage_3ss'})

    df_rmats_const_5ss_temp = df_rmats_const_5ss.reset_index()
    rmats_const_5ss = pbt.BedTool.from_dataframe(df_rmats_const_5ss_temp[['chr', 'start', 'end', 'index', 14, 'strand']])
    rmats_const_5ss = rmats_const_5ss.slop(l=window, r=window, g=fai, s=True).intersect(pbt.BedTool(xl_bed), s=True, wao=True)
    df_const_5ss = rmats_const_5ss.to_dataframe(header=None)
    df_const_5ss = df_const_5ss[[0, 1, 2, 3, 4, 10, 12]]
    df_const_5ss_grouped10 = df_const_5ss.groupby(by=[0, 1, 2, 3, 4])[10].sum().to_frame()
    df_const_5ss_grouped_12 = df_const_5ss.groupby(by=[0, 1, 2, 3, 4])[12].sum().to_frame()
    df_const_5ss = df_const_5ss_grouped10.join(df_const_5ss_grouped_12).reset_index()
    df_const_5ss = df_const_5ss.set_index(3)
    df_rmats_const_5ss_temp = pd.merge(df_rmats_const_5ss_temp, df_const_5ss[[10, 12]], left_on='index', right_index=True)
    del df_rmats_const_5ss_temp['index']
    df_rmats_const_5ss_temp = df_rmats_const_5ss_temp.set_index(14)
    df_rmats_const_5ss_temp = df_rmats_const_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_const = df_rmats_const_5ss_temp.join(df_rmats_const_3ss_temp, lsuffix='5ss', rsuffix='3ss')
    
    df_rmats_const2 = df_rmats_const[['chr5ss', 'start5ss', 'end5ss', 'strand5ss', 'cDNA_coverage_5ss', 'base_coverage_5ss',
        'chr3ss', 'start3ss', 'end3ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]

    df_rmats_const2 = df_rmats_const2.reset_index()

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & (2*window > (df_rmats_const2['start5ss'] - df_rmats_const2['start3ss'])) & \
                    ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) > window),
                    'start'] = df_rmats_const2['start5ss'] - window
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & (2*window > (df_rmats_const2['start5ss'] - df_rmats_const2['start3ss'])) & \
                    ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) > window),
                    'end'] = df_rmats_const2['start3ss'] + window

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & (2*window > (df_rmats_const2['start3ss'] - df_rmats_const2['start5ss'])) & \
                    ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) > window),
                    'start'] = df_rmats_const2['start3ss'] - window
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & (2*window > (df_rmats_const2['start3ss'] - df_rmats_const2['start5ss'])) & \
                    ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) > window),
                    'end'] = df_rmats_const2['start5ss'] + window

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) < window),
                    'start'] = df_rmats_const2['start3ss']
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '+') & ((df_rmats_const2['start5ss'] - df_rmats_const2['start3ss']) < window),
                    'end'] = df_rmats_const2['start5ss']

    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) < window),
                    'start'] = df_rmats_const2['start5ss']
    df_rmats_const2.loc[(df_rmats_const2['strand5ss'] == '-') & ((df_rmats_const2['start3ss'] - df_rmats_const2['start5ss']) < window),
                    'end'] = df_rmats_const2['start3ss']

    df_rmats_const2['len'] = df_rmats_const2['end'] - df_rmats_const2['start']

    df_rmats_const2['diff'] = abs(df_rmats_const2['start3ss'] - df_rmats_const2['start5ss'])

    df_rmats_const2 = df_rmats_const2[['chr5ss', 'start', 'end', 14, 'base_coverage_5ss', 'strand5ss']]

    df_rmats_const2 = df_rmats_const2[~(df_rmats_const2['start'].isna() | df_rmats_const2['end'].isna())]

    df_rmats_const2['start'] = df_rmats_const2['start'].astype(int)
    df_rmats_const2['end'] = df_rmats_const2['end'].astype(int)

    rmats_const2 = pbt.BedTool.from_dataframe(df_rmats_const2)
    rmats_const2 = rmats_const2.intersect(pbt.BedTool(xl_bed), s=True, wao=True)

    df_rmats_const2 = rmats_const2.to_dataframe(names=['chrom', 'start', 'end', 'index', 'base_coverage_5ss', 'strand', 'chrom2', 'start2', 'end2', 'useless', 'overlap_score', 'strand2', 'overlap_base'])

    del df_rmats_const2['base_coverage_5ss']
    del df_rmats_const2['chrom2']
    del df_rmats_const2['start2']
    del df_rmats_const2['end2']
    del df_rmats_const2['useless']
    del df_rmats_const2['strand2']

    df_rmats_const2_grouped10 = df_rmats_const2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_score'].sum().to_frame()
    df_rmats_const2_grouped12 = df_rmats_const2.groupby(by=['chrom', 'start', 'end', 'index', 'strand'])['overlap_base'].sum().to_frame()
    df_rmats_const2 = df_rmats_const2_grouped10.join(df_rmats_const2_grouped12).reset_index()
    df_rmats_const2 = df_rmats_const2.set_index('index')

    df_rmats_const_overlap = df_rmats_const.join(df_rmats_const2)

    df_rmats_const_overlap = df_rmats_const_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'overlap_score', 'overlap_base']].fillna(0)

    df_rmats_const_overlap.loc[df_rmats_const_overlap['overlap_score'] == -1, 'overlap_score'] = 0
    df_rmats_const_overlap.loc[df_rmats_const_overlap['cDNA_coverage_5ss'] == -1, 'cDNA_coverage_5ss'] = 0
    df_rmats_const_overlap.loc[df_rmats_const_overlap['cDNA_coverage_3ss'] == -1, 'cDNA_coverage_3ss'] = 0

    df_rmats_const_overlap['cDNA_coverage_5ss'] = df_rmats_const_overlap['cDNA_coverage_5ss'] - df_rmats_const_overlap['overlap_score']
    df_rmats_const_overlap['cDNA_coverage_3ss'] = df_rmats_const_overlap['cDNA_coverage_3ss'] - df_rmats_const_overlap['overlap_score']
    df_rmats_const_overlap['base_coverage_5ss'] = df_rmats_const_overlap['base_coverage_5ss'] - df_rmats_const_overlap['overlap_base']
    df_rmats_const_overlap['base_coverage_3ss'] = df_rmats_const_overlap['base_coverage_3ss'] - df_rmats_const_overlap['overlap_base']

    df_rmats_const = df_rmats_const_overlap[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss']]
    

    df_rmats_enh['class'] = 'enh'
    df_rmats_sil['class'] = 'sil'
    df_rmats_enhrest['class'] = 'enhrest'
    df_rmats_silrest['class'] = 'silrest'
    df_rmats_ctrl['class'] = 'ctrl'
    df_rmats_const['class'] = 'const'

    df_temp = pd.concat([df_rmats_enh, df_rmats_sil, df_rmats_enhrest, df_rmats_silrest, df_rmats_ctrl, df_rmats_const])
    df_temp = df_temp[['cDNA_coverage_5ss', 'base_coverage_5ss', 'cDNA_coverage_3ss', 'base_coverage_3ss', 'class']]

    df_temp.to_csv(f'{output_dir}/{name}_temp_labeled.tsv', sep='\t', index=None)

    df_out = df_rmats.merge(df_temp, how='left', left_index = True, right_index=True)

    df_out.to_csv(f'{output_dir}/{name}_labeled.tsv', sep='\t', index=None)

    
    if de_source == 'rmats':
        df_final_3ss_upstream = df_enh_3ss_upstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2).merge(
            df_sil_3ss_upstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
            df_ctrl_3ss_upstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True)
        df_final_5ss_upstream = df_enh_5ss_upstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2).merge(
            df_sil_5ss_upstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
            df_ctrl_5ss_upstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True)
        df_final_3ss_downstream = df_enh_3ss_downstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2).merge(
            df_sil_3ss_downstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
            df_ctrl_3ss_downstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True)
        df_final_5ss_downstream = df_enh_5ss_downstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2).merge(
            df_sil_5ss_downstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True).merge(
            df_ctrl_5ss_downstream.rolling(smoothing, center=True, win_type="gaussian").mean(std=2), left_index=True, right_index=True)
        df_final_3ss_upstream.to_csv(f'{output_dir}/{name}_final_3ss_upstream.tsv', sep='\t', index=None)
        df_final_5ss_upstream.to_csv(f'{output_dir}/{name}_final_5ss_upstream.tsv', sep='\t', index=None)
        df_final_3ss_downstream.to_csv(f'{output_dir}/{name}_final_3ss_downstream.tsv', sep='\t', index=None)
        df_final_5ss_downstream.to_csv(f'{output_dir}/{name}_final_5ss_downstream.tsv', sep='\t', index=None)
        
    
    # if z_test:
    #     df_radnom_ctrl_3ss = get_random_coverage(df_coverage_3ss, window, n_exons, n_samples).rename(
    #         columns={'mean': 'mean_ctrl', 'std': 'std_ctrl'})
    #     df_radnom_enh_3ss = get_random_coverage(df_coverage_enh_3ss, window, n_exons, n_samples).rename(
    #         columns={'mean': 'mean_enh', 'std': 'std_enh'})
    #     df_radnom_sil_3ss = get_random_coverage(df_coverage_sil_3ss, window, n_exons, n_samples).rename(
    #         columns={'mean': 'mean_sil', 'std': 'std_sil'})
    #     df_z_score_3ss =  df_enh_3ss.merge(df_sil_3ss, left_index=True, right_index=True).merge(
    #         df_radnom_ctrl_3ss, left_index=True, right_index=True).merge(df_ctrl_3ss, left_index=True, right_index=True).merge(
    #     df_radnom_enh_3ss, left_index=True, right_index=True).merge(df_radnom_sil_3ss, left_index=True, right_index=True)
    #     df_radnom_ctrl_5ss = get_random_coverage(df_coverage_5ss, window, n_exons, n_samples ).rename(
    #         columns={'mean': 'mean_ctrl', 'std': 'std_ctrl'})
    #     df_radnom_enh_5ss = get_random_coverage(df_coverage_enh_5ss, window, n_exons, n_samples).rename(
    #         columns={'mean': 'mean_enh', 'std': 'std_enh'})
    #     df_radnom_sil_5ss = get_random_coverage(df_coverage_sil_5ss, window, n_exons, n_samples).rename(
    #         columns={'mean': 'mean_sil', 'std': 'std_sil'})

    #     df_z_score_5ss =  df_enh_5ss.merge(df_sil_5ss, left_index=True, right_index=True).merge(
    #         df_radnom_ctrl_5ss, left_index=True, right_index=True).merge(df_ctrl_5ss, left_index=True, right_index=True).merge(
    #     df_radnom_enh_5ss, left_index=True, right_index=True).merge(df_radnom_sil_5ss, left_index=True, right_index=True)

    #     df_z_score_3ss['z_test_enh'] = (df_z_score_3ss['mean_enh'] - df_z_score_3ss['mean_ctrl']) / (
    #         df_z_score_3ss['std_enh'].pow(2) / n_exons + df_z_score_3ss['std_ctrl'].pow(2) / n_exons)**(1/2)
    #     df_z_score_3ss['z_test_sil'] = (df_z_score_3ss['mean_sil'] - df_z_score_3ss['mean_ctrl']) / (
    #         df_z_score_3ss['std_sil'].pow(2) / n_exons + df_z_score_3ss['std_ctrl'].pow(2) / n_exons)**(1/2)
    #     df_z_score_5ss['z_test_enh'] = (df_z_score_5ss['mean_enh'] - df_z_score_5ss['mean_ctrl']) / (
    #         df_z_score_5ss['std_enh'].pow(2) / n_exons + df_z_score_5ss['std_ctrl'].pow(2) / n_exons)**(1/2)
    #     df_z_score_5ss['z_test_sil'] = (df_z_score_5ss['mean_sil'] - df_z_score_5ss['mean_ctrl']) / (
    #         df_z_score_5ss['std_sil'].pow(2) / n_exons + df_z_score_5ss['std_ctrl'].pow(2) / n_exons)**(1/2)
        
    #     df_z_score_3ss['z_score_enh'] = (df_z_score_3ss.enhanced - df_z_score_3ss['mean_ctrl']) / df_z_score_3ss['std_ctrl']
    #     df_z_score_3ss['z_score_sil'] = (df_z_score_3ss.silenced - df_z_score_3ss['mean_ctrl']) / df_z_score_3ss['std_ctrl']
    #     df_z_score_3ss['z_score_ctrl'] = (df_z_score_3ss.control - df_z_score_3ss['mean_ctrl']) / df_z_score_3ss['std_ctrl']
        
    #     df_z_score_5ss['z_score_enh'] = (df_z_score_5ss.enhanced - df_z_score_5ss['mean_ctrl']) / df_z_score_5ss['std_ctrl']
    #     df_z_score_5ss['z_score_sil'] = (df_z_score_5ss.silenced - df_z_score_5ss['mean_ctrl']) / df_z_score_5ss['std_ctrl']
    #     df_z_score_5ss['z_score_ctrl'] = (df_z_score_5ss.control - df_z_score_5ss['mean_ctrl']) / df_z_score_5ss['std_ctrl']
        
    
    df_final_3ss['relative_position'] = range(-window, window + 1)
    df_final_3ss = df_final_3ss.set_index('relative_position')
    df_final_5ss['relative_position'] = range(-window, window + 1)
    df_final_5ss = df_final_5ss.set_index('relative_position')
    
    if de_source == 'rmats':
        df_final_3ss_upstream['relative_position'] = range(-window, window + 1)
        df_final_3ss_upstream = df_final_3ss_upstream.set_index('relative_position')
        df_final_5ss_upstream['relative_position'] = range(-window, window + 1)
        df_final_5ss_upstream = df_final_5ss_upstream.set_index('relative_position')

        df_final_3ss_downstream['relative_position'] = range(-window, window + 1)
        df_final_3ss_downstream = df_final_3ss_downstream.set_index('relative_position')
        df_final_5ss_downstream['relative_position'] = range(-window, window + 1)
        df_final_5ss_downstream = df_final_5ss_downstream.set_index('relative_position')
    
    df_fisher_3ss['enh_fold_change_smooth'] = df_fisher_3ss['enh_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['sil_fold_change_smooth'] = df_fisher_3ss['sil_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enh_fold_change_smooth'] = df_fisher_5ss['enh_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['sil_fold_change_smooth'] = df_fisher_5ss['sil_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['enhrest_fold_change_smooth'] = df_fisher_3ss['enhrest_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['silrest_fold_change_smooth'] = df_fisher_3ss['silrest_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_3ss['const_fold_change_smooth'] = df_fisher_3ss['const_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['enhrest_fold_change_smooth'] = df_fisher_5ss['enhrest_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['silrest_fold_change_smooth'] = df_fisher_5ss['silrest_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)
    df_fisher_5ss['const_fold_change_smooth'] = df_fisher_5ss['const_fold_change'].rolling(smoothing, center=True, win_type="gaussian").mean(std=2)

    df_fisher_3ss['relative_position'] = range(-window, window + 1)
    df_fisher_3ss = df_fisher_3ss.set_index('relative_position')
    df_fisher_5ss['relative_position'] = range(-window, window + 1)
    df_fisher_5ss = df_fisher_5ss.set_index('relative_position')

    df_fisher_3ss.to_csv(f'{output_dir}/{name}_fisher_3ss.tsv', sep='\t', index=None)
    df_fisher_5ss.to_csv(f'{output_dir}/{name}_fisher_5ss.tsv', sep='\t', index=None)
    
    sns.set(rc={'figure.figsize':(26, 14)})
    sns.set_style("whitegrid")
    colors_fig1 = [colors_dict['enh'], colors_dict['sil'], colors_dict['ctrl'], colors_dict['enhrest'], 
        colors_dict['silrest'], colors_dict['const']]
    new_labels = [f'enhanced {str(len(df_rmats_enh_3ss))} exons', f'silenced {str(len(df_rmats_sil_3ss))} exons',
                  f'control {str(len(df_rmats_ctrl_3ss))} exons', f'enhanced rest {str(len(df_rmats_enhrest_3ss))} exons',
                  f'silenced rest {str(len(df_rmats_silrest_3ss))} exons', f'constitutive {str(len(df_rmats_const_3ss))} exons']
    fig1, axs = plt.subplots(2, 2, sharey='row')
    axs[0, 0] = sns.lineplot(data=df_final_3ss.loc[list(range(int(-window), int(window * 0.2))), :],
        palette=colors_fig1, ax=axs[0, 0], linewidth=linewidth, dashes=[(1, 0), (1, 0), (1, 0), (2, 2), (2, 2), (2, 2)])
    axs[0, 0].set_title("Coverage around 3'SS event")
    axs[0, 0].set_ylabel('Normalised coverage of crosslinks or peaks')
    axs[0, 0].set_xlabel("Position relative to 3'SS")
    axs[0, 0].legend(new_labels)
    axs[0, 1] = sns.lineplot(data=df_final_5ss.loc[list(range(int(-window * 0.2), int(window))), :], 
        palette=colors_fig1, ax=axs[0, 1], linewidth=linewidth, dashes=[(1, 0), (1, 0), (1, 0), (2, 2), (2, 2), (2, 2)])
    axs[0, 1].set_title("Coverage around 5'SS event")
    axs[0, 1].set_xlabel("Position relative to 5'SS")
    axs[0, 1].legend(new_labels)
    
    #fig1_2, axs = plt.subplots(1, 2, sharey='row')
    colors_fig1_2 = [colors_dict['enh'], colors_dict['sil'], colors_dict['enhrest'], 
        colors_dict['silrest'], colors_dict['const']]
    axs[1, 0] = sns.lineplot(data=df_fisher_3ss.loc[list(range(int(-window), int(window * 0.2))), 
                ['enh_smooth', 'sil_smooth', 'enhrest_smooth', 'silrest_smooth', 'const_smooth']],
                 palette = colors_fig1_2,  ax=axs[1, 0], linewidth=linewidth, 
                dashes=[(1, 0), (1, 0), (2, 2), (2, 2), (2, 2)])
    axs[1, 0].set_ylabel("-log(pval) of coverage compared to control using Fisher's test")
    axs[1, 0].set_xlabel("Position relative to 3'SS")
    axs[1, 1] = sns.lineplot(data=df_fisher_5ss.loc[list(range(int(-window * 0.2), int(window))), 
                ['enh_smooth', 'sil_smooth', 'enhrest_smooth', 'silrest_smooth', 'const_smooth']],
                 palette = colors_fig1_2,  ax=axs[1, 1], linewidth=linewidth, dashes=[(1, 0), (1, 0), (2, 2), (2, 2), (2, 2)])
    axs[1, 1].set_xlabel("Position relative to 5'SS")
    plt.tight_layout()
    fig1.savefig(f'{output_dir}/{name}_fig1.pdf')
    
    if de_source == 'rmats':
        sns.set(rc={'figure.figsize':(20, 6)})
        sns.set_style("whitegrid")
        fig2, axs = plt.subplots(1, 2, sharey='row')
        colors_fig2 = [colors_dict['enh'], colors_dict['sil'], colors_dict['ctrl']]
        new_labels2 = [f'enhanced {str(len(df_rmats_enh_3ss))} exons', f'silenced {str(len(df_rmats_sil_3ss))} exons',
                  f'control {str(len(df_rmats_ctrl_3ss))} exons']
        sns.lineplot(data=df_final_3ss_upstream.loc[list(range(int(-window), int(window * 0.2))), :], palette=colors_fig2, ax=axs[0], 
            linewidth=linewidth, dashes=[(1, 0), (1, 0), (1, 0)])
        sns.lineplot(data=df_final_5ss_upstream.loc[list(range(int(-window * 0.2), int(window))), :], palette=colors_fig2, ax=axs[1], 
            linewidth=linewidth, dashes=[(1, 0), (1, 0), (1, 0)])
        axs[0].set_title("Coverage around 3'SS of upstream exon")
        axs[1].set_title("Coverage around 5'SS of upstream exon")
        axs[0].set_ylabel('Normalised coverage of crosslinks or peaks')
        axs[0].set_xlabel("Position relative to 3'SS")
        axs[1].set_xlabel("Position relative to 5'SS")
        axs[0].legend(new_labels2)
        axs[1].legend(new_labels2)
        plt.tight_layout()
        fig2.savefig(f'{output_dir}/{name}_upstream.pdf')

        sns.set(rc={'figure.figsize':(20, 6)})
        sns.set_style("whitegrid")
        fig3, axs = plt.subplots(1, 2, sharey='row')
        sns.lineplot(data=df_final_3ss_downstream.loc[list(range(int(-window), int(window * 0.2))), :], palette=colors_fig2, ax=axs[0], 
            linewidth=linewidth, dashes=[(1, 0), (1, 0), (1, 0)])
        sns.lineplot(data=df_final_5ss_downstream.loc[list(range(int(-window * 0.2), int(window))), :], palette=colors_fig2, ax=axs[1], 
            linewidth=linewidth, dashes=[(1, 0), (1, 0), (1, 0)])
        axs[0].set_title("Coverage around 3'SS of downstream exon")
        axs[1].set_title("Coverage around 5'SS of downstream exon")
        axs[0].set_ylabel('Normalised coverage of crosslinks or peaks')
        axs[0].set_xlabel("Position relative to 3'SS")
        axs[1].set_xlabel("Position relative to 5'SS")
        axs[0].legend(new_labels2)
        axs[1].legend(new_labels2)
        plt.tight_layout()
        fig3.savefig(f'{output_dir}/{name}_downstream.pdf')
    
    max_b = max(max(enh_3ss_no_xl.values()), max(sil_3ss_no_xl.values()), max(ctrl_3ss_no_xl.values()),
               max(enhrest_3ss_no_xl.values()), max(silrest_3ss_no_xl.values()), max(const_3ss_no_xl.values()))
    ctrl_3ss_y, ctrl_3ss_x = np.histogram(np.array(list(ctrl_3ss_no_xl.values())), range=(0, max_b), density=True)
    enh_3ss_y, enh_3ss_x = np.histogram(np.array(list(enh_3ss_no_xl.values())), range=(0, max_b), density=True)
    sil_3ss_y, sil_3ss_x = np.histogram(np.array(list(sil_3ss_no_xl.values())), range=(0, max_b), density=True)
    ctrl_5ss_y, ctrl_5ss_x = np.histogram(np.array(list(ctrl_5ss_no_xl.values())), range=(0, max_b), density=True)
    enh_5ss_y, enh_5ss_x = np.histogram(np.array(list(enh_5ss_no_xl.values())), range=(0, max_b), density=True)
    sil_5ss_y, sil_5ss_x = np.histogram(np.array(list(sil_5ss_no_xl.values())), range=(0, max_b), density=True)
    
    enhrest_3ss_y, enhrest_3ss_x = np.histogram(np.array(list(enhrest_3ss_no_xl.values())), range=(0, max_b), density=True)
    silrest_3ss_y, silrest_3ss_x = np.histogram(np.array(list(silrest_3ss_no_xl.values())), range=(0, max_b), density=True)
    const_3ss_y, const_3ss_x = np.histogram(np.array(list(const_3ss_no_xl.values())), range=(0, max_b), density=True)
    
    enhrest_5ss_y, enhrest_5ss_x = np.histogram(np.array(list(enhrest_5ss_no_xl.values())), range=(0, max_b), density=True)
    silrest_5ss_y, silrest_5ss_x = np.histogram(np.array(list(silrest_5ss_no_xl.values())), range=(0, max_b), density=True)
    const_5ss_y, const_5ss_x = np.histogram(np.array(list(const_5ss_no_xl.values())), range=(0, max_b), density=True)
    
    sns.set(rc={'figure.figsize':(10, 6)})
    sns.set_style("whitegrid")
    fig3, ax = plt.subplots()
    ax.set_xscale('log', basex=10)
    sns.lineplot(x=ctrl_3ss_x[:-1], y=ctrl_3ss_y, color=colors_dict['ctrl'])
    sns.lineplot(x=enh_3ss_x[:-1], y=enh_3ss_y, color=colors_dict['enh'])
    sns.lineplot(x=sil_3ss_x[:-1], y=sil_3ss_y, color=colors_dict['sil'])
    sns.lineplot(x=enhrest_3ss_x[:-1], y=enhrest_3ss_y, color=colors_dict['enhrest'])
    sns.lineplot(x=silrest_3ss_x[:-1], y=silrest_3ss_y, color=colors_dict['silrest'])
    sns.lineplot(x=const_3ss_x[:-1], y=const_3ss_y, color=colors_dict['const'])
    plt.tight_layout()
    fig3.savefig(f'{output_dir}/{name}_3ss_coverage_dstribution.pdf')
    fig4, ax = plt.subplots()
    ax.set_xscale('log', basex=10)
    sns.lineplot(x=ctrl_5ss_x[:-1], y=ctrl_5ss_y, color=colors_dict['ctrl'])
    sns.lineplot(x=enh_5ss_x[:-1], y=enh_5ss_y, color=colors_dict['enh'])
    sns.lineplot(x=sil_5ss_x[:-1], y=sil_5ss_y, color=colors_dict['sil'])
    sns.lineplot(x=enhrest_5ss_x[:-1], y=enhrest_5ss_y, color=colors_dict['enhrest'])
    sns.lineplot(x=silrest_5ss_x[:-1], y=silrest_5ss_y, color=colors_dict['silrest'])
    sns.lineplot(x=const_5ss_x[:-1], y=const_5ss_y, color=colors_dict['const'])
    plt.tight_layout()
    fig4.savefig(f'{output_dir}/{name}_5ss_coverage_dstribution.pdf')

     
    df_temp_enh = pd.DataFrame.from_dict(enh_no_xl, orient='index')
    flat_list_enh = [item for sublist in df_temp_enh.values for item in sublist]
    df_temp_sil = pd.DataFrame.from_dict(sil_no_xl, orient='index')
    flat_list_sil = [item for sublist in df_temp_sil.values for item in sublist]
    df_temp_ctrl = pd.DataFrame.from_dict(ctrl_no_xl, orient='index')
    flat_list_ctrl = [item for sublist in df_temp_ctrl.values for item in sublist]
    
    df_temp_enhrest = pd.DataFrame.from_dict(enhrest_no_xl, orient='index')
    flat_list_enhrest = [item for sublist in df_temp_enhrest.values for item in sublist]
    df_temp_silrest = pd.DataFrame.from_dict(silrest_no_xl, orient='index')
    flat_list_silrest = [item for sublist in df_temp_silrest.values for item in sublist]
    
    df_temp_const = pd.DataFrame.from_dict(const_no_xl, orient='index')
    flat_list_const = [item for sublist in df_temp_const.values for item in sublist]
    
    df_temp_enh = pd.DataFrame(df_rmats_enh_3ss['score'])
    df_temp_enh['no_xl'] = flat_list_enh
    df_temp_enh['class'] = 'enhanced'
    df_temp_sil = pd.DataFrame(df_rmats_sil_3ss['score'])
    df_temp_sil['no_xl'] = flat_list_sil
    df_temp_sil['class'] = 'silenced'
    df_temp_ctrl = pd.DataFrame(df_rmats_ctrl_3ss['score'])
    df_temp_ctrl['no_xl'] = flat_list_ctrl
    df_temp_ctrl['class'] = 'control'
    
    df_temp_enhrest = pd.DataFrame(df_rmats_enhrest_3ss['score'])
    df_temp_enhrest['no_xl'] = flat_list_enhrest
    df_temp_enhrest['class'] = 'enhanced_rest'
    df_temp_silrest = pd.DataFrame(df_rmats_silrest_3ss['score'])
    df_temp_silrest['no_xl'] = flat_list_silrest
    df_temp_silrest['class'] = 'silenced_rest'
    df_temp_const = pd.DataFrame(df_rmats_const_3ss['score'])
    df_temp_const['no_xl'] = flat_list_const
    df_temp_const['class'] = 'const'
    df_temp_foreground = pd.concat([df_temp_enh, df_temp_sil, df_temp_ctrl, df_temp_const])
    df_temp_background = pd.concat([df_temp_enhrest, df_temp_silrest])

    sns.set(rc={'figure.figsize':(20, 6)})
    sns.set_style("whitegrid")
    sns.scatterplot(x=df_temp_background['score'], y=df_temp_background['no_xl'], hue=df_temp_background['class'], 
        palette=[colors_dict['enhrest'], colors_dict['silrest']], alpha=0.8, s=8)
    sns.scatterplot(x=df_temp_foreground['score'], y=df_temp_foreground['no_xl'], hue=df_temp_foreground['class'], 
        palette=[colors_dict['enh'], colors_dict['sil'], colors_dict['ctrl'], colors_dict['const']], s=8)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/{name}_exon_classes.pdf')
    pbt.helpers.cleanup()

 
if __name__=='__main__':
    (
        de_file,
        xl_bed,
        fai,
        output_folder,
        window,
        smoothing, 
        min_ctrl,
        max_ctrl,
        max_inclusion,
        max_fdr,
        max_enh,
        min_sil
    ) = cli()
    
    run_rna_map(de_file, xl_bed, fai, window, smoothing, 
        min_ctrl, max_ctrl, max_inclusion, max_fdr, max_enh, min_sil, output_folder)


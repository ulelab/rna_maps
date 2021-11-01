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

sns.set_style("whitegrid")
colors_dict = {'all': '#D3D3D3', 'ctrl': '#408F76', 'enh': '#F30C08', 'sil': '#005299', 'enhrest': '#FFB122', 'silrest': '#6DC2F5', 'const': '#666666'}
linewidth = 3
dashes = False


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


# def get_random_coverage(df, window, n_exons = 150, n_samples = 300):
#     event_starts = list(df[df.position == 1].index)
#     event_len = 2 * window + 1
#     positions = range(event_len)
#     print(f'{n_samples} coverages for  randomly sampled {n_exons} exons')
#     random_coverages = pd.DataFrame(index= range(2 * window + 1))
#     times = []
#     counter = 0
#     for i in range(n_samples):
#         random_coverages_temp = pd.DataFrame(index= range(2 * window + 1))
#         start_time = time.time()
#         random_events = random.choices(event_starts, k=n_exons)
#         indices = [[i + j for j in positions] for i in random_events]
#         for j, index in enumerate(indices):
#             df_coverage = df.loc[index, :].rename(columns={'coverage': f'coverage{i}_{j}'})
#             if df_coverage[f'coverage{i}_{j}'].sum() == 0:
#                 continue
#             else:
#                 counter += 1
#                 df_plot = df_coverage.copy()
#                 df_plot.loc[df_plot.strand=='+', 'map'] = df_plot['position']
#                 df_plot.loc[df_plot.strand=='-', 'map'] = abs(2 * window + 2  - df_plot['position'])
#                 df_plot.map = df_plot.map.astype('int32')
#                 df_plot = df_plot[[f'coverage{i}_{j}', 'map']]
#                 df_plot = df_plot.set_index('map')
#                 random_coverages_temp = random_coverages_temp.join(df_plot)
#         random_coverages_temp = random_coverages_temp.sum(axis=1) / n_exons
#         random_coverages = random_coverages.join(pd.DataFrame(random_coverages_temp).rename(columns= {0: f'coverage_{i}'}))
#         end_time = time.time()
#         times.append(end_time - start_time)
#     random_coverages['mean'] = random_coverages.mean(axis=1)
#     random_coverages['std'] = random_coverages.std(axis=1)
#     random_coverages['new_index'] = random_coverages.index + 1
#     random_coverages.set_index('new_index', inplace=True)
#     return random_coverages[['mean', 'std']]


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

def get_3ss5ss_exons_whippet(df_exons):
    exons_stranded = df_exons.copy()
    exons_stranded.loc[df_exons['strand'] == '-', 'start'] = df_exons['end'] # replace start with end
    exons_stranded.loc[df_exons['strand'] == '-', 'end'] = df_exons['start'] # replace end with start
    exons_3ss = exons_stranded[['chrom', 'start', 'DeltaPsi', 'Probability', 'strand', 'exon_len', 'inclusion', '14']].copy() # takes all the starts for 3ss
    exons_5ss = exons_stranded[['chrom', 'end', 'DeltaPsi', 'Probability', 'strand', 'exon_len', 'inclusion', '14']].copy() # takes all the ends for 5ss
    exons_3ss['end'] = exons_3ss['start'] + 1 # calculates end from start
    exons_5ss['start'] = exons_5ss['end'] # assigns start to end for 5ss
    exons_5ss['end'] = exons_5ss['start'] + 1 # calulates end from start
    return exons_3ss[['chrom', 'start', 'end', 'DeltaPsi', 'Probability', 'strand', 'exon_len', 'inclusion', '14']], exons_5ss[
        ['chrom', 'start', 'end', 'DeltaPsi', 'Probability', 'strand', 'exon_len', 'inclusion', '14']]


def run_rna_map(de_file, xl_bed, fai, window=300, smoothing=15, 
        min_ctrl=-0.05, max_ctrl=0.05, max_inclusion=0.99, max_fdr=0.1, max_enh=-0.05, min_sil=0.05, min_prob_whippet=0.9, output_dir='.',
       #n_exons = 150, n_samples = 300, z_test=False
       ):
    name = de_file.split('/')[-1].replace('.txt', '').replace('.gz', '')
    df_fai = pd.read_csv(fai, sep='\t', header=None)
    chroms = set(df_fai[0].values)
    rmats = pd.read_csv(de_file, sep='\t')
    if 'exonStart_0base' in rmats.columns:
        rmats = rmats[rmats['chr'].isin(chroms)]
        rmats['inclusion'] = (rmats.IncLevel1.str.split(',') + rmats.IncLevel2.str.split(','))
        rmats['inclusion'] = rmats['inclusion'].apply(lambda x: sum([float(y) for y in x if y != 'NA']) / len(x))
        df_rmats =  rmats.loc[ : ,['chr', 'exonStart_0base', 'exonEnd', 'FDR', 'IncLevelDifference', 'strand', 'inclusion', 
                                   'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']].rename(
            columns={'chr': 0, 'exonStart_0base': 1, 'exonEnd': 2, 'FDR': 3, 'IncLevelDifference': 4, 'strand': 5, 'inclusion': 6, 
                     'upstreamES': 7, 'upstreamEE': 8, 'downstreamES': 9, 'downstreamEE': 10})
        print('Using rmats output file')
        de_source = 'rmats'
        df_rmats[11] = df_rmats[2] - df_rmats[1] # calculate exon length
        df_rmats[12] = df_rmats[8] - df_rmats[7] # calculate exon length
        df_rmats[13] = df_rmats[10] - df_rmats[9] #  calculate exon length
        
    elif 'DeltaPsi' in rmats.columns:
        rmats = pd.read_csv(de_file, sep='\t', index_col=False)
        df_rmats = rmats[rmats['Type'] == 'CE']
        df_rmats['chrom'] = df_rmats['Coord'].apply(lambda x: x.split(':')[0])
        df_rmats = df_rmats[df_rmats['chrom'].isin(chroms)]
        df_rmats['start'] = df_rmats['Coord'].apply(lambda x: x.split(':')[1].split('-')[0])
        df_rmats['start'] = df_rmats['start'].astype(int)
        df_rmats['end'] = df_rmats['Coord'].apply(lambda x: x.split(':')[1].split('-')[1])
        df_rmats['end'] = df_rmats['end'].astype(int)
        df_rmats = df_rmats.rename(columns={'Strand': 'strand'})
        print('Using whippet output file')
        de_source = 'whippet'
        df_rmats['inclusion'] = (df_rmats['Psi_A'] + df_rmats['Psi_B']) / 2
        df_rmats['exon_len'] = df_rmats['end'] - df_rmats['start'] # calculate exon length
    else:
        df_rmats = pd.read_csv(de_file, sep='\t', header=None, index_col=None, 
                               dtype={0: 'str', 1: 'int', 2: 'int', 3: 'str', 4: 'float', 5: 'str'})
        try:
            df_rmats[3] = df_rmats[3].astype(float) # check if there is a numerical value that can be used for filtering e.g. FDR
        except: # if not assume all events are to be used
            del df_rmats[3]
            df_rmats[3] = 0.0 #non rmats files will probably not have FDR in 4th column so this value will make pass all
            df_rmats[6] = 0.0 #non rmats files will probably not have inclusion in 7th column so this value will make pass all
            df_rmats = df_rmats[[0, 1, 2, 3, 4, 5, 6]]

    df_rmats = df_rmats.reset_index()
    df_rmats = df_rmats.rename(columns={'index': 14})

    if de_source == 'rmats':
        rmats_3ss, rmats_5ss = get_3ss5ss_exons(df_rmats)
        print(f'There are: {len(rmats_3ss)} 3ss and {len(rmats_5ss)} 5ss')
        col_rename = {0: 'chr', 1: 'start', 2: 'end', 3: 'name', 4: 'score', 5: 'strand', 6:'inclusion', 7: 'upstream_es', 
            8: 'upstream_ee', 9: 'downstream_es', 10: 'downstream_ee', 11: 'exon_len', 12: 'upstream_exon_len', 13: 'downstream_exon_len'}
        df_rmats_ctrl_3ss = rmats_3ss[(rmats_3ss[4] > min_ctrl) & (rmats_3ss[4] < max_ctrl) &
                                      (rmats_3ss[6] < max_inclusion) & (rmats_3ss[3] > max_fdr)].rename(columns=col_rename)
        df_rmats_ctrl_5ss = rmats_5ss[(rmats_5ss[4] > min_ctrl) & (rmats_5ss[4] < max_ctrl) &
                                      (rmats_5ss[6] < max_inclusion) & (rmats_5ss[3] > max_fdr)].rename(columns=col_rename)
        df_rmats_const_3ss = rmats_3ss[(rmats_3ss[4] > min_ctrl) & (rmats_3ss[4] < max_ctrl) &
                                      (rmats_3ss[6] > max_inclusion) & (rmats_3ss[3] > max_fdr)].rename(columns=col_rename)
        df_rmats_const_5ss = rmats_5ss[(rmats_5ss[4] > min_ctrl) & (rmats_5ss[4] < max_ctrl) &
                                      (rmats_5ss[6] > max_inclusion) & (rmats_5ss[3] > max_fdr)].rename(columns=col_rename)
        del df_rmats_ctrl_3ss['inclusion']; del df_rmats_ctrl_5ss['inclusion']; del rmats_3ss[6]; del rmats_5ss[6]
        df_rmats_enh_3ss = rmats_3ss[(rmats_3ss[4] < max_enh) & (rmats_3ss[3] < max_fdr)].rename(columns=col_rename)
        df_rmats_enh_5ss = rmats_5ss[(rmats_5ss[4] < max_enh) & (rmats_5ss[3] < max_fdr)].rename(columns=col_rename) 
        df_rmats_sil_3ss = rmats_3ss[(rmats_3ss[4] > min_sil) & (rmats_3ss[3] < max_fdr)].rename(columns=col_rename)
        df_rmats_sil_5ss = rmats_5ss[(rmats_5ss[4] > min_sil) & (rmats_3ss[3] < max_fdr)].rename(columns=col_rename)
        df_rmats_enhrest_3ss = rmats_3ss[(rmats_3ss[4] < max_enh) & (rmats_3ss[3] >= max_fdr)].rename(columns=col_rename)
        df_rmats_enhrest_5ss = rmats_5ss[(rmats_5ss[4] < max_enh) & (rmats_5ss[3] >= max_fdr)].rename(columns=col_rename) 
        df_rmats_silrest_3ss = rmats_3ss[(rmats_3ss[4] > min_sil) & (rmats_3ss[3] >= max_fdr)].rename(columns=col_rename)
        df_rmats_silrest_5ss = rmats_5ss[(rmats_5ss[4] > min_sil) & (rmats_3ss[3] >= max_fdr)].rename(columns=col_rename)

        df_rmats_ctrl_3ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_ctrl_5ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_const_3ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_const_5ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_enh_3ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_enh_5ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_sil_3ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_sil_5ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_enhrest_3ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_enhrest_5ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_silrest_3ss.drop_duplicates(subset=['start', 'end'], inplace=True)
        df_rmats_silrest_5ss.drop_duplicates(subset=['start', 'end'], inplace=True)

        df_rmats_sil_3ss = df_rmats_sil_3ss[~(
            ((df_rmats_sil_3ss['start'].isin(df_rmats_enh_3ss['start'])) & (df_rmats_sil_3ss['end'].isin(df_rmats_enh_3ss['end'])))
            )]
        df_rmats_sil_5ss = df_rmats_sil_5ss[~(
            ((df_rmats_sil_5ss['start'].isin(df_rmats_enh_5ss['start'])) & (df_rmats_sil_5ss['end'].isin(df_rmats_enh_5ss['end'])))
            )]
        
        df_rmats_enh_3ss = df_rmats_enh_3ss[~(
            ((df_rmats_enh_3ss['start'].isin(df_rmats_sil_3ss['start'])) & (df_rmats_enh_3ss['end'].isin(df_rmats_sil_3ss['end'])))
            )]
        df_rmats_enh_5ss = df_rmats_enh_5ss[~(
            ((df_rmats_enh_5ss['start'].isin(df_rmats_sil_5ss['start'])) & (df_rmats_enh_5ss['end'].isin(df_rmats_sil_5ss['end'])))
            )]

        df_rmats_silrest_3ss = df_rmats_silrest_3ss[~(
            ((df_rmats_silrest_3ss['start'].isin(df_rmats_enhrest_3ss['start'])) & (df_rmats_silrest_3ss['end'].isin(df_rmats_enhrest_3ss['end']))) |\
            ((df_rmats_silrest_3ss['start'].isin(df_rmats_sil_3ss['start'])) & (df_rmats_silrest_3ss['end'].isin(df_rmats_sil_3ss['end']))) |\
            ((df_rmats_silrest_3ss['start'].isin(df_rmats_enh_3ss['start'])) & (df_rmats_silrest_3ss['end'].isin(df_rmats_enh_3ss['end'])))
            )]
        df_rmats_silrest_5ss = df_rmats_silrest_5ss[~(
            ((df_rmats_silrest_5ss['start'].isin(df_rmats_enhrest_5ss['start'])) & (df_rmats_silrest_5ss['end'].isin(df_rmats_enhrest_5ss['end']))) |\
            ((df_rmats_silrest_5ss['start'].isin(df_rmats_sil_5ss['start'])) & (df_rmats_silrest_5ss['end'].isin(df_rmats_sil_5ss['end']))) |\
            ((df_rmats_silrest_5ss['start'].isin(df_rmats_enh_5ss['start'])) & (df_rmats_silrest_5ss['end'].isin(df_rmats_enh_5ss['end'])))
            )]
        
        df_rmats_enhrest_3ss = df_rmats_enhrest_3ss[~(
            ((df_rmats_enhrest_3ss['start'].isin(df_rmats_enhrest_3ss['start'])) & (df_rmats_enhrest_3ss['end'].isin(df_rmats_silrest_3ss['end']))) |\
            ((df_rmats_enhrest_3ss['start'].isin(df_rmats_sil_3ss['start'])) & (df_rmats_enhrest_3ss['end'].isin(df_rmats_sil_3ss['end']))) |\
            ((df_rmats_enhrest_3ss['start'].isin(df_rmats_enh_3ss['start'])) & (df_rmats_enhrest_3ss['end'].isin(df_rmats_enh_3ss['end'])))
            )]
        df_rmats_enhrest_5ss = df_rmats_enhrest_5ss[~(
            ((df_rmats_enhrest_5ss['start'].isin(df_rmats_enhrest_5ss['start'])) & (df_rmats_enhrest_5ss['end'].isin(df_rmats_silrest_5ss['end']))) |\
            ((df_rmats_enhrest_5ss['start'].isin(df_rmats_sil_5ss['start'])) & (df_rmats_enhrest_5ss['end'].isin(df_rmats_sil_5ss['end']))) |\
            ((df_rmats_enhrest_5ss['start'].isin(df_rmats_enh_5ss['start'])) & (df_rmats_enhrest_5ss['end'].isin(df_rmats_enh_5ss['end'])))
            )]

        df_rmats_ctrl_3ss = df_rmats_ctrl_3ss[~(
            ((df_rmats_ctrl_3ss['start'].isin(df_rmats_sil_3ss['start'])) & (df_rmats_ctrl_3ss['end'].isin(df_rmats_sil_3ss['end']))) |\
            ((df_rmats_ctrl_3ss['start'].isin(df_rmats_enh_3ss['start'])) & (df_rmats_ctrl_3ss['end'].isin(df_rmats_enh_3ss['end']))) |\
            ((df_rmats_ctrl_3ss['start'].isin(df_rmats_silrest_3ss['start'])) & (df_rmats_ctrl_3ss['end'].isin(df_rmats_silrest_3ss['end']))) |\
            ((df_rmats_ctrl_3ss['start'].isin(df_rmats_enhrest_3ss['start'])) & (df_rmats_ctrl_3ss['end'].isin(df_rmats_enhrest_3ss['end'])))
            )]
        df_rmats_ctrl_5ss = df_rmats_ctrl_5ss[~(
            ((df_rmats_ctrl_5ss['start'].isin(df_rmats_sil_5ss['start'])) & (df_rmats_ctrl_5ss['end'].isin(df_rmats_sil_5ss['end']))) |\
            ((df_rmats_ctrl_5ss['start'].isin(df_rmats_enh_5ss['start'])) & (df_rmats_ctrl_5ss['end'].isin(df_rmats_enh_5ss['end']))) |\
            ((df_rmats_ctrl_5ss['start'].isin(df_rmats_silrest_5ss['start'])) & (df_rmats_ctrl_5ss['end'].isin(df_rmats_silrest_5ss['end']))) |\
            ((df_rmats_ctrl_5ss['start'].isin(df_rmats_enhrest_5ss['start'])) & (df_rmats_ctrl_5ss['end'].isin(df_rmats_enhrest_5ss['end'])))
            )]

        df_rmats_const_3ss = df_rmats_const_3ss[~(
            ((df_rmats_const_3ss['start'].isin(df_rmats_sil_3ss['start'])) & (df_rmats_const_3ss['end'].isin(df_rmats_sil_3ss['end']))) |\
            ((df_rmats_const_3ss['start'].isin(df_rmats_enh_3ss['start'])) & (df_rmats_const_3ss['end'].isin(df_rmats_enh_3ss['end']))) |\
            ((df_rmats_const_3ss['start'].isin(df_rmats_silrest_3ss['start'])) & (df_rmats_const_3ss['end'].isin(df_rmats_silrest_3ss['end']))) |\
            ((df_rmats_const_3ss['start'].isin(df_rmats_enhrest_3ss['start'])) & (df_rmats_const_3ss['end'].isin(df_rmats_enhrest_3ss['end']))) |\
            ((df_rmats_const_3ss['start'].isin(df_rmats_ctrl_3ss['start'])) & (df_rmats_const_3ss['end'].isin(df_rmats_ctrl_3ss['end'])))
            )]
        df_rmats_const_5ss = df_rmats_const_5ss[~(
            ((df_rmats_const_5ss['start'].isin(df_rmats_sil_5ss['start'])) & (df_rmats_const_5ss['end'].isin(df_rmats_sil_5ss['end']))) |\
            ((df_rmats_const_5ss['start'].isin(df_rmats_enh_5ss['start'])) & (df_rmats_const_5ss['end'].isin(df_rmats_enh_5ss['end']))) |\
            ((df_rmats_const_5ss['start'].isin(df_rmats_silrest_5ss['start'])) & (df_rmats_const_5ss['end'].isin(df_rmats_silrest_5ss['end']))) |\
            ((df_rmats_const_5ss['start'].isin(df_rmats_enhrest_5ss['start'])) & (df_rmats_const_5ss['end'].isin(df_rmats_enhrest_5ss['end']))) |\
            ((df_rmats_const_5ss['start'].isin(df_rmats_ctrl_5ss['start'])) & (df_rmats_const_5ss['end'].isin(df_rmats_ctrl_5ss['end'])))
            )]

        index_selected = df_rmats_ctrl_3ss.index.union(df_rmats_enh_3ss.index.union(
            df_rmats_sil_3ss.index.union(df_rmats_enhrest_3ss.index.union(df_rmats_silrest_3ss.index.union(df_rmats_const_3ss.index)))))

    elif de_source == 'whippet':
        rmats_3ss, rmats_5ss = get_3ss5ss_exons_whippet(df_rmats)
        print(f'There are: {len(rmats_3ss)} 3ss and {len(rmats_5ss)} 5ss')
        df_rmats_ctrl_3ss = rmats_3ss[(rmats_3ss['DeltaPsi'] > min_ctrl) & (rmats_3ss['DeltaPsi'] < max_ctrl) &
                                      (rmats_3ss['inclusion'] < max_inclusion) & (rmats_3ss['Probability'] > max_fdr)]
        df_rmats_ctrl_5ss = rmats_5ss[(rmats_5ss['DeltaPsi'] > min_ctrl) & (rmats_5ss['DeltaPsi'] < max_ctrl) &
                                      (rmats_5ss['inclusion'] < max_inclusion) & (rmats_5ss['Probability'] > max_fdr)]
        df_rmats_const_3ss = rmats_3ss[(rmats_3ss['DeltaPsi'] > min_ctrl) & (rmats_3ss['DeltaPsi'] < max_ctrl) &
                                      (rmats_3ss['inclusion'] > max_inclusion) & (rmats_3ss['Probability'] > max_fdr)]
        df_rmats_const_5ss = rmats_5ss[(rmats_5ss['DeltaPsi'] > min_ctrl) & (rmats_5ss['DeltaPsi'] < max_ctrl) &
                                      (rmats_5ss['inclusion'] > max_inclusion) & (rmats_5ss['Probability'] > max_fdr)]
        del df_rmats_ctrl_3ss['inclusion']; del df_rmats_ctrl_5ss['inclusion']; del rmats_3ss['inclusion']; del rmats_5ss['inclusion']
        df_rmats_enh_3ss = rmats_3ss[(rmats_3ss['DeltaPsi'] < max_enh) & (rmats_3ss['Probability'] > min_prob_whippet)]
        df_rmats_enh_5ss = rmats_5ss[(rmats_5ss['DeltaPsi'] < max_enh) & (rmats_5ss['Probability'] > min_prob_whippet)] 
        df_rmats_sil_3ss = rmats_3ss[(rmats_3ss['DeltaPsi'] > min_sil) & (rmats_3ss['Probability'] > min_prob_whippet)]
        df_rmats_sil_5ss = rmats_5ss[(rmats_5ss['DeltaPsi'] > min_sil) & (rmats_3ss['Probability'] > min_prob_whippet)]
        df_rmats_enhrest_3ss = rmats_3ss[(rmats_3ss['DeltaPsi'] < max_enh) & (rmats_3ss['Probability'] <= min_prob_whippet)]
        df_rmats_enhrest_5ss = rmats_5ss[(rmats_5ss['DeltaPsi'] < max_enh) & (rmats_5ss['Probability'] <= min_prob_whippet)] 
        df_rmats_silrest_3ss = rmats_3ss[(rmats_3ss['DeltaPsi'] > min_sil) & (rmats_3ss['Probability'] <= min_prob_whippet)]
        df_rmats_silrest_5ss = rmats_5ss[(rmats_5ss['DeltaPsi'] > min_sil) & (rmats_3ss['Probability'] <= min_prob_whippet)]


        index_selected = df_rmats_ctrl_3ss.index.union(df_rmats_enh_3ss.index.union(
            df_rmats_sil_3ss.index.union(df_rmats_enhrest_3ss.index.union(df_rmats_silrest_3ss.index.union(df_rmats_const_3ss.index)))))
                
#     print(f'There are {len(df_rmats_enh_3ss)} enhanced exons, {len(df_rmats_sil_3ss)} silenced exons, ' \
#           f'{len(df_rmats_ctrl_3ss)} control exons and {len(df_rmats_rest_3ss)} left over exons')
    if len(df_rmats_ctrl_3ss) == 0:
        print('No control exons, try changing filtering parameters or file')
        return
    if len(df_rmats_enh_3ss) == 0 and len(df_rmats_sil_3ss) == 0:
        print('No regulated exons, try changing filtering parameters or file')
        return
    if de_source == 'rmats':
        exon_len_all = df_rmats[11]
    elif de_source == 'whippet':
        exon_len_all = df_rmats['exon_len']
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
        upstream_exon_len_all = df_rmats[12]
        upstream_exon_len_ctrl = df_rmats_ctrl_3ss['upstream_exon_len']
        upstream_exon_len_enh = df_rmats_enh_3ss['upstream_exon_len']
        upstream_exon_len_sil = df_rmats_sil_3ss['upstream_exon_len']
        upstream_exon_len_enhrest = df_rmats_enhrest_3ss['upstream_exon_len']
        upstream_exon_len_silrest = df_rmats_silrest_3ss['upstream_exon_len']
        upstream_exon_len_const = df_rmats_const_3ss['upstream_exon_len']
        df_upstream_exon_len = pd.DataFrame(data={'enhanced': upstream_exon_len_enh, 'enhanced_rest': upstream_exon_len_enhrest, 'silenced': upstream_exon_len_sil, 
            'silenced_rest': upstream_exon_len_silrest, 'control': upstream_exon_len_ctrl, 
            'constitutive': upstream_exon_len_const, 'all_exons': upstream_exon_len_all})

        downstream_exon_len_all = df_rmats[13]
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
    if de_source == 'whippet':
        rename_cols = {'chrom': 'chr', 'DeltaPsi': 'name', 'Probability': 'score'}
        df_rmats_enh_3ss = df_rmats_enh_3ss.rename(columns=rename_cols)
        df_rmats_enh_5ss = df_rmats_enh_5ss.rename(columns=rename_cols)
        df_rmats_sil_3ss = df_rmats_sil_3ss.rename(columns=rename_cols)
        df_rmats_sil_5ss = df_rmats_sil_5ss.rename(columns=rename_cols)
        df_rmats_ctrl_3ss = df_rmats_ctrl_3ss.rename(columns=rename_cols)
        df_rmats_ctrl_5ss = df_rmats_ctrl_5ss.rename(columns=rename_cols)
        df_rmats_enhrest_3ss = df_rmats_enhrest_3ss.rename(columns=rename_cols)
        df_rmats_enhrest_5ss = df_rmats_enhrest_5ss.rename(columns=rename_cols)
        df_rmats_silrest_3ss = df_rmats_silrest_3ss.rename(columns=rename_cols)
        df_rmats_silrest_5ss = df_rmats_silrest_5ss.rename(columns=rename_cols)
        df_rmats_const_3ss = df_rmats_const_3ss.rename(columns=rename_cols)
        df_rmats_const_5ss = df_rmats_const_5ss.rename(columns=rename_cols)
    
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
    df_rmats_enh_3ss_temp = df_rmats_enh_3ss_temp[[10, 12]]
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
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss_temp[[10, 12]]
    df_rmats_enh_5ss_temp = df_rmats_enh_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_enh = df_rmats_enh_5ss_temp.join(df_rmats_enh_3ss_temp)

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
    df_rmats_sil_3ss_temp = df_rmats_sil_3ss_temp[[10, 12]]
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

    df_rmats_sil = df_rmats_sil_5ss_temp.join(df_rmats_sil_3ss_temp)

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
    df_rmats_enhrest_3ss_temp = df_rmats_enhrest_3ss_temp[[10, 12]]
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
    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss_temp[[10, 12]]
    df_rmats_enhrest_5ss_temp = df_rmats_enhrest_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_enhrest = df_rmats_enhrest_5ss_temp.join(df_rmats_enhrest_3ss_temp)

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
    df_rmats_silrest_3ss_temp = df_rmats_silrest_3ss_temp[[10, 12]]
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
    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss_temp[[10, 12]]
    df_rmats_silrest_5ss_temp = df_rmats_silrest_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_silrest = df_rmats_silrest_5ss_temp.join(df_rmats_silrest_3ss_temp)

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
    df_rmats_ctrl_3ss_temp = df_rmats_ctrl_3ss_temp[[10, 12]]
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
    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss_temp[[10, 12]]
    df_rmats_ctrl_5ss_temp = df_rmats_ctrl_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_ctrl = df_rmats_ctrl_5ss_temp.join(df_rmats_ctrl_3ss_temp)

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
    df_rmats_const_3ss_temp = df_rmats_const_3ss_temp[[10, 12]]
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
    df_rmats_const_5ss_temp = df_rmats_const_5ss_temp[[10, 12]]
    df_rmats_const_5ss_temp = df_rmats_const_5ss_temp.rename(columns={10: 'cDNA_coverage_5ss', 12: 'base_coverage_5ss'})

    df_rmats_const = df_rmats_const_5ss_temp.join(df_rmats_const_3ss_temp)

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

    # df_fisher_3ss.to_csv(f'{output_dir}/{name}_fisher_3ss.tsv', sep='\t', index=None)
    # df_fisher_5ss.to_csv(f'{output_dir}/{name}_fisher_5ss.tsv', sep='\t', index=None)
    
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
    import sys

    de_file = sys.argv[1]
    xl_bed = sys.argv[2]
    fai = sys.argv[3]
    output_folder= sys.argv[4]
    window = int(sys.argv[5])#300
    smoothing = int(sys.argv[6])#15 
    min_ctrl = float(sys.argv[7])#-0.05
    max_ctrl = float(sys.argv[8])#0.05
    max_inclusion= float(sys.argv[9])#0.9
    max_fdr = float(sys.argv[10])#0.1
    max_enh = float(sys.argv[11])#-0.05
    min_sil = float(sys.argv[12])#0.05
    min_prob_whippet = float(sys.argv[13])#0.9
    
    run_rna_map(de_file, xl_bed, fai, window, smoothing, 
        min_ctrl, max_ctrl, max_inclusion, max_fdr, max_enh, min_sil, min_prob_whippet, output_folder)

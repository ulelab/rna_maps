import pandas as pd

def load_and_filter_rmats(de_file, fai):

    # 1.1 Extract valid chromosomes from .fai
    df_fai = pd.read_csv(fai, sep="\t", header=None, usecols=[0])
    valid_chroms = set(df_fai[0].tolist())

    # 1.2 Read rMATS file
    rmats_df = pd.read_csv(de_file, sep="\t")
    # 1.3 Filter by chromosome
    if 'chr' not in rmats_df.columns:
        raise ValueError("Expected 'chr' column in rMATS file.")
    rmats_df = rmats_df[rmats_df['chr'].isin(valid_chroms)].copy()
    rmats_df.reset_index(drop=True, inplace=True)
    return rmats_df

def compute_max_psi_and_dedup(rmats_df):
    df = rmats_df.copy()

    # 2.1 Compute maxPSI: merge IncLevel1 and IncLevel2 strings, split by comma, ignore 'NA'
    def calc_max_psi(row):
        psi_vals = []
        for col in ['IncLevel1', 'IncLevel2']:
            if pd.isna(row[col]):
                continue

            psi_vals.extend([float(x) for x in row[col].split(',') if x != 'NA'])
            # print(col,psi_vals)
        return max(psi_vals) if psi_vals else 0.0

    df['maxPSI'] = df.apply(calc_max_psi, axis=1)

    # 2.2 Rename IncLevelDifference to dPSI if not already
    if 'IncLevelDifference' in df.columns:
        df.rename(columns={'IncLevelDifference': 'dPSI'}, inplace=True)
    if 'dPSI' not in df.columns:
        raise ValueError("Expected 'IncLevelDifference' or 'dPSI' column in rMATS file.")

    # 2.3 Select relevant columns for downstream analysis
    df_subset = df[['chr', 'exonStart_0base', 'exonEnd', 'strand',
                    'FDR', 'dPSI', 'maxPSI',
                    'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']].copy()
    df_subset.reset_index(drop=True, inplace=True)

    # 2.4 De-duplicate by (chr, exonStart_0base, exonEnd, strand), keeping the row with highest |dPSI|
    df_subset['abs_dPSI'] = df_subset['dPSI'].abs()
    df_subset['rank'] = df_subset.groupby(
        ['chr', 'exonStart_0base', 'exonEnd', 'strand']
    )['abs_dPSI'] .transform(lambda x: abs(x).rank(ascending=False)) < 2
    df_clean = df_subset[df_subset['rank'] == 1].drop(columns=['abs_dPSI', 'rank']).reset_index(drop=True)
    return df_clean

def get_ss_bed(df, pos_col, neg_col):
    df['exon_id'] = df['category'] + "_" + df['chr'].astype(str) + ":" + df['exonStart_0base'].astype(str) + "-" + df['exonEnd'].astype(str) + ";" + df['strand'].astype(str)

    ss_pos = df.loc[df['strand'] == "+", ['chr', pos_col, pos_col, 'exon_id', 'FDR', 'strand']]
    ss_pos.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_pos.start = ss_pos.start.transform(lambda x: x-1)

    ss_n = df.loc[df['strand'] == "-", ['chr', neg_col, neg_col, 'exon_id', 'FDR', 'strand']]
    ss_n.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_n.end = ss_n.end.transform(lambda x: x+1)

    ss = pd.concat([ss_pos, ss_n])

    return ss

def prepare_ss_bed(df):
    ss_map = {
        'upstream_3ss': get_ss_bed(df, 'upstreamES', 'upstreamEE'),
        'upstream_5ss': get_ss_bed(df, 'upstreamEE', 'upstreamES'),
        'middle_3ss':   get_ss_bed(df, 'exonStart_0base', 'exonEnd'),
        'middle_5ss':   get_ss_bed(df, 'exonEnd', 'exonStart_0base'),
        'downstream_3ss': get_ss_bed(df, 'downstreamES', 'downstreamEE'),
        'downstream_5ss': get_ss_bed(df, 'downstreamEE', 'downstreamES')
    }
    return ss_map

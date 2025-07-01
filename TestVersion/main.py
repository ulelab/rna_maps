from utils.data import load_and_filter_rmats, compute_max_psi_and_dedup, prepare_ss_bed
from utils.analysis import categorize_exons, downsample_categories, compute_exon_lengths, smooth_coverage, get_coverage_plot, do_multivalency_analysis
from utils.visualize import plot_exon_length_boxplot, plot_ss_heatmap, plot_enrichment_curves
import pandas as pd
import os



multivalency = False
de_file= 'test/chr21_PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt'
fai = 'test/homosapien-hg37-chr21.fa.fai'
genome_fasta = 'test/homosapien-hg37-chr21.fa'
xl_bed = 'test/chr21_hela_ptbp1_iclip_sorted_merged.bed'
output_dir = 'test_output'
os.makedirs(output_dir,exist_ok=True)
germsdir = os.getcwd()
window = 300
# colors_dict = {'all': '#D3D3D3', 'ctrl': '#408F76', 'enh': '#F30C08', 'sil': '#005299', 'enhrest': '#FFB122', 'silrest': '#6DC2F5', 'const': '#666666'}
# smoothing = 15


def main():
    rmats = load_and_filter_rmats(de_file, fai)
    df_rmats = compute_max_psi_and_dedup(rmats)
    df_cat = categorize_exons(df_rmats)
    # print(df_cat)
    df_sampled, exon_counts_after, exon_counts_before = downsample_categories(df_cat,0)
    # print(df_sampled)
    # print(exon_counts_after)
    # print(exon_counts_before)

    length_df = compute_exon_lengths(df_sampled)
    plot_exon_length_boxplot(length_df,
                             output_path=f"{output_dir}/exon_length.pdf",
                             title="Exon Length Distribution",
                             titles = ["Upstream Exon", "Middle Exon", "Downstream exon"])
    # print(length_df)
    ss_map = prepare_ss_bed(df_sampled)
    # print(ss_map)


    exon_counts_after = pd.Series(exon_counts_after)
    line_dict = {}
    heat_dict = {}
    if xl_bed is not None:
        for label, bed_df in ss_map.items():
            line_df, heat_df = get_coverage_plot(
                xl_bed, bed_df, fai, window, exon_counts_after, label
            )
            # 二值化并平滑 heat_df
            # heat_df['coverage'] = (heat_df['coverage'] > 0).astype(int)
            heat_df = smooth_coverage(heat_df)

            line_dict[label] = line_df
            heat_dict[label] = heat_df

        plotting_df = pd.concat([line_dict[l] for l in ss_map.keys()], ignore_index=True)
        plotting_df.to_csv(f'{output_dir}/RNAmap.tsv', sep="\t", index=False)


        combined_heat = pd.concat([heat_dict[l] for l in ss_map.keys()], ignore_index=True)
        exon_totals = combined_heat.groupby('exon_id')['coverage'].sum().reset_index().set_index('exon_id')
        exon_names = combined_heat[['exon_id','name']].drop_duplicates().set_index('exon_id')['name']

        exon_info = exon_totals.copy()
        exon_info['name'] = exon_names
        exon_info = exon_info.sort_values(['name', 'coverage'], ascending=[True, False])
        sorted_ids = exon_info.index.tolist()


        plot_ss_heatmap(
            label_data_dict=heat_dict,
            exon_totals=exon_totals,
            exon_names=exon_names,
            sorted_exon_ids=sorted_ids,
            labels=list(ss_map.keys()),
            window=window,
            output_path=f'{output_dir}/heatmap.pdf',
            title="Splice-site Coverage Heatmap"
        )


        plot_enrichment_curves(
            plotting_df=plotting_df,
            labels=list(ss_map.keys()),
            window=window,
            exon_categories_after=exon_counts_after,
            exon_categories_before=exon_counts_before,
            output_path=f"{output_dir}/RNAmap_-log10pvalue.pdf",
            title="Splice-site -log10(p-value) Enrichment"
        )

        if multivalency:
            do_multivalency_analysis(
                ss_map=ss_map,
                fai=fai,
                window=window,
                genome_fasta=genome_fasta,
                output_dir=output_dir,
                FILEname="RNAmap",
                germsdir=germsdir
            )

if __name__=='__main__':
    main()
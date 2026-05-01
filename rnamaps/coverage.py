"""CLIP coverage calculation around splice sites with Fisher's exact test."""

import numpy as np
import pandas as pd
import pybedtools as pbt
import scipy.stats as stats


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

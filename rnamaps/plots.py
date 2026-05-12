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


def _smooth_rows_gaussian(matrix, window_size=10, std=2):
    """Vectorised per-row Gaussian rolling-mean smoothing along columns.

    Mirrors the previous per-(exon, label) Python loop but runs in C
    via pandas' ``DataFrame.rolling``. Edge NaNs are filled with the
    original (unsmoothed) value to match the prior behaviour.
    """
    if matrix.size == 0:
        return matrix.astype(np.float32, copy=False)
    df = pd.DataFrame(matrix.T.astype(np.float32, copy=False))
    smoothed = df.rolling(
        window=window_size, center=True, win_type='gaussian'
    ).mean(std=std)
    raw = np.array(df.values, copy=True)
    sm = np.array(smoothed.values, copy=True)
    mask = np.isnan(sm)
    sm[mask] = raw[mask]
    return sm.T


def plot_heatmap(region_covs, exon_categories, window, all_sites,
                 output_dir, FILEname):
    """Generate per-exon heatmap from per-region coverage matrices.

    Parameters
    ----------
    region_covs : list of rnamaps.coverage.RegionCoverage
        One per splice-site region. Each carries a per-exon × per-position
        coverage matrix.
    """
    # Total exons covered table: count distinct exons with any coverage
    # within each (label, category) pair.
    count_rows = []
    for rc in region_covs:
        if rc.n_exons == 0:
            continue
        any_hit = (rc.matrix > 0).any(axis=1)
        if not any_hit.any():
            continue
        names_with_hit = rc.exon_names[any_hit]
        unique, counts = np.unique(names_with_hit, return_counts=True)
        for cat, cnt in zip(unique, counts):
            count_rows.append({'label': rc.label,
                               'name': cat,
                               'exon_count': int(cnt)})
    count_heat_df = pd.DataFrame(
        count_rows, columns=['label', 'name', 'exon_count']
    )
    final_heat_df = count_heat_df.pivot(
        index='name', columns='label', values='exon_count'
    ).fillna(0).reset_index() if not count_heat_df.empty else pd.DataFrame(
        {'name': []}
    )
    exon_categories_df = exon_categories.reset_index()
    exon_categories_df.columns = ['name', 'total_exons_after_subsetting']
    final_heat_df = final_heat_df.merge(
        exon_categories_df, on='name', how='left'
    )
    final_heat_df.to_csv(
        f'{output_dir}/{FILEname}_totalExonsCovered.tsv', sep="\t", index=False
    )

    if not all_sites:
        labels = ['upstream_5ss', 'middle_3ss', 'middle_5ss', 'downstream_3ss']
    else:
        labels = ['upstream_3ss', 'upstream_5ss', 'middle_3ss', 'middle_5ss',
                  'downstream_3ss', 'downstream_5ss']

    # Build binarised + smoothed matrices keyed by label, plus a global
    # per-exon "total signal" used to rank rows. All operations stay in
    # numpy so we never materialise the (n_exons × n_positions × n_regions)
    # long-form DataFrame the old path needed.
    cov_by_label = {rc.label: rc for rc in region_covs}
    if not cov_by_label:
        logging.info("No regions for heatmap — skipping")
        return

    # Use the first available region to learn the global exon order /
    # identity. All regions are built from the same parent exon frame so
    # exon_ids align row-wise across regions.
    any_rc = region_covs[0]
    exon_ids = any_rc.exon_ids
    exon_names = any_rc.exon_names
    n_exons = exon_ids.size
    if n_exons == 0:
        logging.info("No exons with signal for heatmap — skipping")
        return

    total_signal = np.zeros(n_exons, dtype=np.float64)
    smoothed_by_label = {}
    for label in labels:
        if label not in cov_by_label:
            continue
        rc = cov_by_label[label]
        bin_mat = (rc.matrix > 0).astype(np.float32, copy=False)
        sm = _smooth_rows_gaussian(bin_mat, window_size=10, std=2)
        smoothed_by_label[label] = sm
        total_signal += sm.sum(axis=1)

    keep = total_signal > 0
    if not keep.any():
        logging.info("No exons with signal for heatmap — skipping")
        return

    keep_idx = np.flatnonzero(keep)
    kept_names = exon_names[keep_idx]
    kept_totals = total_signal[keep_idx]

    # Sort by category, then descending total signal within each category.
    order = np.lexsort((-kept_totals, kept_names))
    sorted_idx = keep_idx[order]
    sorted_exon_ids = exon_ids[sorted_idx].tolist()

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

    names = kept_names[order]
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

    # Name group labels (contiguous runs since rows are sorted by name).
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
        if label not in smoothed_by_label:
            ax = fig.add_subplot(gs[0, i + 1])
            ax.set_facecolor('none')
            ax.text(0.5, 0.5, f"No data for {label}", ha='center', va='center')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(label)
            continue

        rc = cov_by_label[label]
        smoothed = smoothed_by_label[label]
        positions = rc.positions

        if '3ss' in label:
            min_pos, max_pos = 0, window + 50
        else:
            min_pos, max_pos = window - 50, window * 2

        col_mask = (positions >= min_pos) & (positions <= max_pos)
        display_matrix = smoothed[np.ix_(sorted_idx, col_mask)]

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


def _layout(all_sites):
    if not all_sites:
        col_order = ["upstream_5ss", "middle_3ss",
                     "middle_5ss", "downstream_3ss"]
        titles = ["Upstream 5'SS", "Middle 3'SS",
                  "Middle 5'SS", "Downstream 3'SS"]
        col_wrap = 4
    else:
        col_order = ["upstream_3ss", "upstream_5ss",
                     "middle_3ss", "middle_5ss",
                     "downstream_3ss", "downstream_5ss"]
        titles = ["Upstream 3'SS", "Upstream 5'SS",
                  "Middle 3'SS", "Middle 5'SS",
                  "Downstream 3'SS", "Downstream 5'SS"]
        col_wrap = 6
    return col_order, titles, col_wrap


_DEFAULT_HUE_ORDER = [
    'constitutive', 'control', 'enhanced', 'silenced',
    'enhanced_rest', 'silenced_rest',
]


def _decorate_exon_intron(g, col_order, window):
    """Draw exon/intron rectangles, set xlims and tick labels per panel."""
    rect_fraction = 1 / ((window + 50) / 50)
    last_rect = None
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
        last_rect = rect
    return last_rect


def _hue_palette(plotting_df):
    hue_order = [c for c in _DEFAULT_HUE_ORDER
                 if c in plotting_df['name'].unique()]
    palette_map = {
        'constitutive': colors_dict['const'],
        'control': colors_dict['ctrl'],
        'enhanced': colors_dict['enh'],
        'silenced': colors_dict['sil'],
        'enhanced_rest': colors_dict['enhrest'],
        'silenced_rest': colors_dict['silrest'],
    }
    palette = {h: palette_map.get(h, '#999999') for h in hue_order}
    return hue_order, palette


def plot_rna_map(plotting_df, exon_categories, original_counts,
                 window, all_sites, output_dir, FILEname,
                 pvalue_method='fisher', n_perm=None, y_axis='log10p',
                 plot_kind='line', y_col=None, ci_cols=None,
                 clusters_df=None, method_name=None, ylabel=None,
                 subtitle=None):
    """Generate RNA map plots.

    Parameters
    ----------
    plot_kind : {'line', 'ribbon', 'clusters'}
        - ``'line'``: single y-column line plot (default; legacy behaviour).
        - ``'ribbon'``: line plot with a shaded CI band from
          ``ci_cols=(lo_col, hi_col)``.
        - ``'clusters'``: line plot for ``y_col`` plus horizontal markers
          beneath the x-axis at significant clusters from ``clusters_df``.
    y_col : str, optional
        Column to plot on the y-axis. For ``plot_kind='line'`` it
        defaults from ``y_axis`` (legacy behaviour).
    ci_cols : tuple of (str, str), optional
        Lower / upper CI column names for ``plot_kind='ribbon'``.
    clusters_df : DataFrame, optional
        Cluster summary table for ``plot_kind='clusters'``.
    method_name : str, optional
        Output PDF suffix; e.g. ``'log2fc'`` ⇒ ``..._RNAmap_log2fc.pdf``.
        Defaults to ``'-log10pvalue'`` (legacy filename).
    ylabel : str, optional
        Y-axis label override.
    subtitle : str, optional
        Small grey subtitle drawn above the figure.
    """
    sns.set(rc={'figure.figsize': (7, 5)})
    sns.set_style("whitegrid")

    col_order, titles, col_wrap = _layout(all_sites)

    if y_col is None:
        y_col = 'zscore' if y_axis == 'zscore' else '-log10pvalue_smoothed'
    if method_name is None:
        method_name = '-log10pvalue'

    hue_order, palette = _hue_palette(plotting_df)

    if plot_kind == 'line':
        g = sns.relplot(
            data=plotting_df, x='position', y=y_col,
            hue='name', col='label', facet_kws={"sharex": False},
            kind='line', col_wrap=col_wrap, height=5, aspect=4 / 5,
            col_order=col_order, hue_order=hue_order, palette=palette,
        )
    elif plot_kind in ('ribbon', 'clusters'):
        g = sns.FacetGrid(
            plotting_df, col='label', col_order=col_order,
            col_wrap=col_wrap, height=5, aspect=4 / 5,
            hue='name', hue_order=hue_order, palette=palette,
            sharex=False, sharey=True,
        )
        if plot_kind == 'ribbon':
            if ci_cols is None:
                ci_cols = (f"{y_col}_lo", f"{y_col}_hi")
            lo_col, hi_col = ci_cols

            def _ribbon(data, color=None, **kwargs):
                d = data.sort_values('position')
                plt.fill_between(
                    d['position'].values,
                    d[lo_col].values, d[hi_col].values,
                    color=color, alpha=0.2, linewidth=0,
                )
                plt.plot(d['position'].values, d[y_col].values,
                         color=color, linewidth=2)

            g.map_dataframe(_ribbon)
        else:
            def _line(data, color=None, **kwargs):
                d = data.sort_values('position')
                plt.plot(d['position'].values, d[y_col].values,
                         color=color, linewidth=2)
            g.map_dataframe(_line)
        g.add_legend()
    else:
        raise ValueError(
            f"Unknown plot_kind={plot_kind!r}; expected one of "
            f"'line', 'ribbon', 'clusters'."
        )

    for ax, title in zip(g.axes.flat, titles):
        ax.set_title(title)
        ax.axhline(y=0, color='k', alpha=0.2, linewidth=0.5)
        fig = plt.gcf()
        marker_ax = add_enrichment_marker(fig, ax)

    g.set(xlabel='')

    if ylabel is None:
        if y_axis == 'zscore':
            ylabel = 'signed permutation z-score vs control'
        elif pvalue_method == 'permutation':
            ylabel = 'signed -log10(empirical p) vs control'
        else:
            ylabel = '-log10(p value) enrichment / control'
    g.axes[0].set_ylabel(ylabel)

    sns.move_legend(
        g, "upper right",
        bbox_to_anchor=(1, 2),
        ncol=1, title=None, frameon=False,
    )
    leg = g._legend
    set_legend_text(leg, exon_categories, original_counts)

    last_rect = _decorate_exon_intron(g, col_order, window)

    # Cluster bars under each panel.
    if plot_kind == 'clusters' and clusters_df is not None and not clusters_df.empty:
        for i, ss_type in enumerate(col_order):
            ax = g.axes[i]
            sub = clusters_df[
                (clusters_df['label'] == ss_type)
                & (clusters_df['cluster_pvalue'] <= 0.05)
            ]
            if sub.empty:
                continue
            ymin, ymax = ax.get_ylim()
            yrange = ymax - ymin
            for j, (_, row) in enumerate(sub.iterrows()):
                cat = row['name']
                color = palette.get(cat, '#999999')
                yoff = ymin - 0.04 * yrange - 0.02 * yrange * j
                ax.plot(
                    [row['start_pos'], row['end_pos']],
                    [yoff, yoff],
                    color=color, linewidth=4, solid_capstyle='butt',
                    clip_on=False,
                )

    plt.subplots_adjust(wspace=0.05)

    auto_subtitle = None
    if subtitle is None and pvalue_method == 'permutation' and n_perm is not None:
        auto_subtitle = f"Label-permutation test (B={n_perm})."
    final_subtitle = subtitle if subtitle is not None else auto_subtitle
    if final_subtitle:
        g.fig.suptitle(final_subtitle, y=1.02, fontsize=8, color='dimgray')

    out_path = f'{output_dir}/{FILEname}_RNAmap_{method_name}.pdf'
    plt.savefig(
        out_path,
        bbox_extra_artists=([leg, last_rect, marker_ax]),
        bbox_inches='tight', pad_inches=0.8,
    )
    logging.info(f"Saved RNA map to {out_path}")
    plt.close('all')
    pbt.helpers.cleanup()


def plot_roc_curves(roc_curves_df, region_auc_df, all_sites,
                    output_dir, FILEname,
                    aggregators_to_plot=('mean',),
                    xl_score_mode=None):
    """Render per-region ROC curves (one PDF per aggregator).

    Parameters
    ----------
    roc_curves_df : pd.DataFrame
        Columns ``name, label, aggregator, fpr, tpr``. One row per
        (region, category, aggregator, threshold-step).
    region_auc_df : pd.DataFrame or None
        Headline AUC table (one row per region × category × aggregator).
        Used to print AUC in the legend.
    all_sites : bool
        4- vs 6-panel layout, matching ``plot_rna_map``.
    aggregators_to_plot : iterable of str
        Subset of ``{'mean', 'max'}`` to render. Each gets its own PDF.
    xl_score_mode : str, optional
        BED-score mode that produced these ROC curves. When provided
        (multi-mode sweep), the output PDF filenames embed it as
        ``..._xlscore-{mode}_curves_{aggregator}.pdf``; single-mode
        runs keep the legacy ``..._curves_{aggregator}.pdf`` names.
    """
    if roc_curves_df is None or roc_curves_df.empty:
        logging.info("No ROC curves to plot — skipping")
        return

    col_order, titles, col_wrap = _layout(all_sites)
    auc_lookup = {}
    if region_auc_df is not None and not region_auc_df.empty:
        for _, row in region_auc_df.iterrows():
            key = (row['label'], row['name'], row['aggregator'])
            auc_lookup[key] = float(row['auc'])

    for aggregator in aggregators_to_plot:
        sub = roc_curves_df[roc_curves_df['aggregator'] == aggregator]
        if sub.empty:
            continue

        sns.set(rc={'figure.figsize': (10, 8)})
        sns.set_style("whitegrid")
        _, palette = _hue_palette(sub)

        n_panels = len(col_order)
        n_rows = (n_panels + col_wrap - 1) // col_wrap
        fig, axes = plt.subplots(
            n_rows, col_wrap,
            figsize=(4 * col_wrap, 4 * n_rows),
            squeeze=False,
        )
        axes_flat = axes.flatten()

        for i, region_label in enumerate(col_order):
            ax = axes_flat[i]
            region_sub = sub[sub['label'] == region_label]
            ax.plot([0, 1], [0, 1], color='lightgrey',
                    linestyle='--', linewidth=1)
            if region_sub.empty:
                ax.set_title(f"{titles[i]} (no data)")
            else:
                for cat, cat_df in region_sub.groupby('name'):
                    cat_df = cat_df.sort_values('fpr')
                    auc_val = auc_lookup.get(
                        (region_label, cat, aggregator), None
                    )
                    label_txt = (f"{cat} (AUC={auc_val:.3f})"
                                 if auc_val is not None
                                 else cat)
                    ax.plot(
                        cat_df['fpr'].values, cat_df['tpr'].values,
                        color=palette.get(cat, '#999999'),
                        linewidth=2, label=label_txt,
                    )
                ax.legend(loc='lower right', fontsize=8, frameon=False)
                ax.set_title(titles[i])
            ax.set_xlim([-0.01, 1.01])
            ax.set_ylim([-0.01, 1.01])
            ax.set_xlabel('False positive rate')
            ax.set_ylabel('True positive rate')

        for j in range(n_panels, len(axes_flat)):
            axes_flat[j].set_visible(False)

        fig.suptitle(
            f"ROC curves (window-{aggregator} per-exon score)",
            fontsize=10, color='dimgray',
        )
        plt.tight_layout(rect=[0, 0, 1, 0.97])

        mode_suffix = f"_xlscore-{xl_score_mode}" if xl_score_mode else ""
        out_path = (
            f'{output_dir}/{FILEname}_RNAmap_roc_auc'
            f'{mode_suffix}_curves_{aggregator}.pdf'
        )
        plt.savefig(out_path, bbox_inches='tight')
        logging.info(f"Saved ROC curves to {out_path}")
        plt.close(fig)

import matplotlib
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

colors_dict = {'all': '#D3D3D3', 'ctrl': '#408F76', 'enh': '#F30C08', 'sil': '#005299', 'enhrest': '#FFB122', 'silrest': '#6DC2F5', 'const': '#666666'}

def plot_exon_length_boxplot(
    length_df,
    output_path,
    title = None,
    titles = None
):

    sns.set_style("whitegrid")
    # sns.set(rc={'figure.figsize': (15, 5)})
    # plt.figure(figsize=(15, 5))
    col_categories = ["upstream_exon_length", "regulated_exon_length", "downstream_exon_length"]
    order_categories = ["control", "constitutive", "enhanced", "enhanced_rest", "silenced", "silenced_rest"]
    palette = [colors_dict['ctrl'], colors_dict['const'], colors_dict['enh'],colors_dict['enhrest'], colors_dict['sil'], colors_dict['silrest'], colors_dict['all']]
    # palette = sns.color_palette("Set2", n_colors=len(order_categories))

    g = sns.catplot(
        data=length_df,
        x='category',
        y='exon_length',
        col='exon_type',
        kind='box',
        col_wrap=3,
        showfliers=False,
        col_order = col_categories,
        order=order_categories,
        palette=palette,
        hue='category',
        legend=False
    )
    if title:
        g.fig.suptitle(title, y=1, fontsize=10)
    if titles:
        for ax, subtitle in zip(g.axes.flat, titles):
            ax.set_title(subtitle)

    g.set_axis_labels("", "Exon Length (bp)")
    for ax in g.axes.flatten():
        ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    # plt.show()
    g.savefig(output_path, dpi=300)
    plt.close()

def plot_ss_heatmap(
        label_data_dict,
        exon_totals,
        exon_names: pd.Series,
        sorted_exon_ids,
        labels,
        window,
        output_path,
        title=None
):
    num_exons = len(sorted_exon_ids)
    unique_names = exon_names.loc[sorted_exon_ids].unique().tolist()
    name_to_idx = {nm: i for i, nm in enumerate(unique_names)}

    # Build name_matrix
    name_matrix = np.zeros((num_exons, 1), dtype=int)
    for i, eid in enumerate(sorted_exon_ids):
        name_matrix[i, 0] = name_to_idx[exon_names.loc[eid]]

    # Figure and GridSpec
    num_labels = len(labels)
    width = max(15, num_labels * 4)
    height = max(3, num_exons * 0.002)
    fig = plt.figure(figsize=(width, height))
    gs = GridSpec(1, num_labels + 1, width_ratios=[1] + [3] * num_labels, figure=fig)

    # 8.1 Plot name column heatmap
    ax_names = fig.add_subplot(gs[0, 0])
    # name_cmap = plt.get_cmap('tab20', len(unique_names))
    color_palette = plt.cm.tab10.colors[:len(unique_names)]
    name_cmap = LinearSegmentedColormap.from_list('name_cmap', [(1, 1, 1)] + list(color_palette),
                                                  N=len(unique_names) + 1)
    sns.heatmap(
        name_matrix,
        ax=ax_names,
        cmap=name_cmap,
        cbar=False,
        xticklabels=False,
        yticklabels=False,
        linewidths=0,
        rasterized=True
    )
    # Annotate groups and draw horizontal lines between groups
    start_idx = 0
    prev_name = exon_names.loc[sorted_exon_ids[0]]
    for i, eid in enumerate(sorted_exon_ids + [None]):
        if i == num_exons or exon_names.loc[eid] != prev_name:
            end_idx = i - 1
            mid = (start_idx + end_idx) / 2
            ax_names.text(
                x=0.5, y=mid, s=prev_name, va='center', ha='center',
                rotation=0, fontsize=10, fontweight='bold'
            )
            if i < num_exons:
                ax_names.hlines(i, *ax_names.get_xlim(), colors='white', linewidth=0.5)
            if i < num_exons:
                prev_name = exon_names.loc[eid]
                start_idx = i

    ax_names.set_title("Exon Groups", fontsize=10)
    ax_names.axis('off')

    # 8.2 Plot heatmap for each splice-site label
    for col_idx, label in enumerate(labels):
        ax = fig.add_subplot(gs[0, col_idx + 1])
        if label not in label_data_dict or label_data_dict[label].empty:
            ax.text(0.5, 0.5, f"No data for\n{label}", ha='center', va='center', fontsize=10)
            ax.axis('off')
            continue

        df_label = label_data_dict[label]
        # Pivot to get rows=exon_id, cols=position, values=coverage (0/1)
        pivot = df_label.pivot(index='exon_id', columns='position', values='coverage').fillna(0)
        # Determine x-axis positions to display
        if '3ss' in label:
            pos_min, pos_max = 0, window + 50
        else:  # '5ss'
            pos_min, pos_max = window - 50, window * 2
        sel_positions = [p for p in pivot.columns if pos_min <= p <= pos_max]
        # Build display matrix: shape (num_exons, len(sel_positions))
        # display_matrix = np.zeros((num_exons, len(sel_positions)), dtype=int)
        # for i, eid in enumerate(sorted_exon_ids):
        #    if eid in pivot.index:
        #       display_matrix[i, :] = pivot.loc[eid, sel_positions].values
        display_matrix = pivot.loc[sorted_exon_ids, sel_positions].fillna(0).values

        vmax_val = df_label['coverage'].max()

        sns.heatmap(
            display_matrix,
            ax=ax,
            cmap=colormaps['viridis'],
            vmin=0,
            vmax=vmax_val,
            cbar=False,
            xticklabels=False,
            yticklabels=False,
            linewidths=0,
            rasterized=True
        )
        # Draw horizontal lines between groups
        for i, eid in enumerate(sorted_exon_ids[:-1]):
            if exon_names.loc[eid] != exon_names.loc[sorted_exon_ids[i + 1]]:
                ax.hlines(i + 1, *ax.get_xlim(), colors='white', linewidth=0.5)

        ax.set_title(label.replace('_', ' ').title(), fontsize=10)
        ax.axis('off')

    if title:
        fig.suptitle(title, fontsize=12, y=1.02)
    plt.subplots_adjust(left=0.15)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_enrichment_curves(
        plotting_df,
        labels,
        window,
        exon_categories_after,
        exon_categories_before,
        output_path,
        title=None
):
    sns.set_style("whitegrid")
    g = sns.relplot(
        data=plotting_df,
        x='position',
        y='-log10pvalue_smoothed',
        hue='name',
        col='label',
        kind='line',
        col_wrap=len(labels),
        height=4,
        aspect=1.2,
        facet_kws={"sharex": False},
        col_order=labels
    )

    # Customize each subplot
    rect_fraction = 1 / ((window + 50) / 50)
    for i, (ax, label) in enumerate(zip(g.axes.flat, labels)):
        # 1)
        if i % 2 == 0:  # 3'SS panels
            ax.set_xlim(0, window - 50)
            ticks = np.arange(0, window + 50 + 1, 50)
        else:  # 5'SS panels
            ax.set_xlim(window - 50, 2 * window)
            ticks = np.arange(window - 50, 2 * window + 1, 50)

        labels_txt = [
            "" if t in (ticks[0], ticks[-1]) else str(int(t - window))
            for t in ticks
        ]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels_txt)

        # 2) axhline
        ax.axhline(0, color='lightgray', linewidth=0.5)

        # 3) set_title
        ax.set_title(label.replace('_', ' ').title(), fontsize=10)
        ax.set_xlabel('')
        ax.set_ylabel('-log10(p-value)')

        # 4) rect1
        color1 = 'midnightblue' if (i == 2 or i == 3) else 'slategrey'
        rect1 = matplotlib.patches.Rectangle(
            (1 - rect_fraction, -0.2) if i % 2 == 0 else (0, -0.2),
            width=rect_fraction,
            height=0.1,
            color=color1, alpha=1,
            transform=ax.transAxes, clip_on=False
        )
        ax.add_patch(rect1)

        #    rect2
        rect2 = matplotlib.patches.Rectangle(
            (0, -0.15) if i % 2 == 0 else (rect_fraction, -0.15),
            width=1 - rect_fraction,
            height=0.001,
            color='slategrey', alpha=1,
            transform=ax.transAxes, clip_on=False
        )
        ax.add_patch(rect2)

    # Adjust legend: show exon group names and counts
    legend = g._legend
    legend.set_title("Exon Groups")
    # Build custom legend text
    legend_texts = []
    for cat, cnt in exon_categories_after.items():
        orig_cnt = exon_categories_before.get(cat, 0)
        legend_texts.append(f"{cat}({cnt}, subset from {orig_cnt})")
    # Place text inside legend title area
    # legend.texts[0].set_text("\n".join(legend_texts))
    # legend.set_title('\n'.join(legend_texts))
    for i, legend_text in enumerate(legend_texts):
        legend.texts[i].set_text(legend_text)
    sns.move_legend(
        g,
        loc='upper left',
        bbox_to_anchor=(1.02, 0.8),
        ncol=1,
        frameon=False
    )
    plt.subplots_adjust(right=0.75)

    if title:
        plt.subplots_adjust(top=0.9)
        g.fig.suptitle(title, fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

"""Shared plotting helpers: legend formatting and enrichment direction marker."""

import pandas as pd


def set_legend_text(legend, exon_categories, original_counts=None):
    """Set legend text with optional subset information."""
    exon_cat = pd.DataFrame({
        'name': exon_categories.index,
        'number_exons': exon_categories.values
    })
    categories = exon_cat['name'].unique()
    legend.set_bbox_to_anchor([1.08, 0.75])
    legend.set_title("")

    legend_idx = 0
    for category in ['constitutive', 'control', 'enhanced', 'silenced']:
        if category in categories:
            count = exon_cat[exon_cat['name'] == category]['number_exons'].values[0]
            text = f"{category.capitalize()} ({count}"
            if original_counts and category in original_counts:
                orig_count = original_counts[category]
                if orig_count > count:
                    text += f", subset from {orig_count}"
            text += ")"
            legend.texts[legend_idx].set_text(text)
            legend_idx += 1


def add_enrichment_marker(fig, ax):
    """Add enrichment/depletion marker with arrows."""
    y_min, y_max = ax.get_ylim()
    marker_height = 0.4
    normalized_bottom = (0 - y_min) / (y_max - y_min) - marker_height / 2

    marker_ax = fig.add_axes(
        [0.94, normalized_bottom, 0.06, marker_height],
        label=f"enrichment_marker_{ax.get_title()}"
    )
    marker_ax.set_xticks([])
    marker_ax.set_yticks([])
    for spine in marker_ax.spines.values():
        spine.set_visible(False)

    marker_ax.text(
        0.5, 0.5, "↑\nenriched\n\n\n↓\ndepleted",
        ha='center', va='center', fontsize=11,
        transform=marker_ax.transAxes
    )
    return marker_ax

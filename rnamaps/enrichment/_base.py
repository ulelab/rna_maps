"""Shared types and helpers for enrichment methods."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Tuple

import pandas as pd


@dataclass
class EnrichmentResult:
    """Per-region output from an enrichment method.

    Attributes
    ----------
    plot_df : pd.DataFrame
        Per (category, position) rows for one splice-site region. Always
        carries ``name``, ``position``, ``label`` plus method-specific
        score columns named in ``y_columns``.
    clusters_df : pd.DataFrame
        Optional cluster-level output (e.g. for cluster-based permutation).
        Empty for methods that don't produce clusters.
    plot_kind : str
        Hint to the plot layer: one of ``"line"``, ``"ribbon"``,
        ``"clusters"``, ``"roc_auc"``.
    y_columns : tuple of str
        Columns in ``plot_df`` to plot. ``"line"`` plots a single y-column;
        ``"ribbon"`` plots ``y_columns`` and looks for matching ``_lo`` /
        ``_hi`` companion columns; ``"clusters"`` plots a single y-column
        plus cluster bars from ``clusters_df``.
    method_name : str
        Short identifier used to name the output PDF and TSV files.
    ylabel : str
        Y-axis label for the plot.
    extras : dict
        Optional method-specific side-data that doesn't fit the standard
        per-(category, position) schema (e.g. ROC curve coordinates for
        ``roc_auc``). The pipeline aggregates this across regions when
        concatenating per-region results.
    """

    plot_df: pd.DataFrame
    clusters_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    plot_kind: str = "line"
    y_columns: Tuple[str, ...] = ("-log10pvalue_smoothed",)
    method_name: str = "enrichment"
    ylabel: str = "enrichment"
    extras: Dict[str, Any] = field(default_factory=dict)

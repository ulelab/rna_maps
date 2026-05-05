"""Enrichment methods for RNA maps.

Each method exposes a ``compute(df_per_exon, exon_categories, label, *, rng,
**opts) -> EnrichmentResult`` function. The pipeline calls one or more
methods (selected via ``--enrichment``) and dispatches their results to
the plot layer based on ``EnrichmentResult.plot_kind``.

Available methods:

- ``bootstrap_contrast`` (default) — per-position delta and log2 fold change
  with bootstrap CIs. Recommended for shape comparison and class-imbalanced
  designs (small regulated, large control).
- ``cluster_perm`` — Maris-Oostenveld cluster-mass permutation. Significance
  per peak with FWER control across positions.
- ``permutation_z`` — legacy label-permutation z-score (in ``rnamaps.permutation``).
- ``fisher`` — legacy per-position Fisher's exact test (in ``rnamaps.coverage``).
"""

from rnamaps.enrichment._base import EnrichmentResult

__all__ = ["EnrichmentResult"]

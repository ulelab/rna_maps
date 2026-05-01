#!/usr/bin/env python3
"""
RNA Maps - entry script.

Generates RNA maps showing RBP binding patterns around regulated exons.
The implementation lives in the `rnamaps` package; this script is the
command line entry point and is invoked as `python rna_maps.py ...`.

Supports two input modes:
  1. rMATS mode:  Takes rMATS differential splicing output, auto-categorises exons
  2. VastDB mode: Takes pre-curated VastDB EVENT ID lists + EVENT_INFO annotation

Usage examples:

  # rMATS mode
  python rna_maps.py -i rMATS_SE.txt -x CLIP.bed -f hg19.fa -fi hg19.fa.fai -o output -p PTBP1

  # VastDB mode
  python rna_maps.py --vastdb_mode \
    --vastdb_enhanced enhanced_ids.txt \
    --vastdb_silenced silenced_ids.txt \
    --vastdb_control control_ids.txt \
    --vastdb_constitutive constitutive_ids.txt \
    --vastdb_annotation EVENT_INFO-hg38.tab \
    -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p AQR_K562
"""

# Apply matplotlib/seaborn config side effects before any plotting code runs.
import rnamaps.config  # noqa: F401
from rnamaps.cli import cli
from rnamaps.pipeline import run_rna_map


if __name__ == '__main__':
    args = cli()
    run_rna_map(args)

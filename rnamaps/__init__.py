"""RNA Maps package.

Generates RNA maps showing RBP binding patterns around regulated exons.
See the top-level rna_maps.py script for the CLI entry point.
"""

from rnamaps.cli import cli
from rnamaps.pipeline import run_rna_map

__all__ = ["cli", "run_rna_map"]

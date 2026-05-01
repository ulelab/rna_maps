"""RNA Maps package.

Generates RNA maps showing RBP binding patterns around regulated exons.
Install with ``pip install -e .`` and run via ``rnamaps`` or ``python -m rnamaps``.
"""

from rnamaps.cli import cli
from rnamaps.pipeline import run_rna_map

__all__ = ["cli", "run_rna_map"]

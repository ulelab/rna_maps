"""Allow running the package with ``python -m rnamaps`` or the ``rnamaps`` command."""

from rnamaps.cli import cli
from rnamaps.pipeline import run_rna_map


def main():
    args = cli()
    run_rna_map(args)


if __name__ == "__main__":
    main()

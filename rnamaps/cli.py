"""Argparse-based command line interface."""

import argparse
import os


def cli():
        optional.add_argument(
            '--y_axis', type=str, default='log10p', choices=['log10p', 'zscore'],
            help="Y-axis for RNA map plots: 'log10p' (default) for signed -log10(p), 'zscore' for signed permutation z-score.")
    parser = argparse.ArgumentParser(
        prog='rnamaps',
        description='Plot CLIP crosslinks around regulated exons to study '
                    'position-dependent impact on pre-mRNA splicing.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input modes:

  rMATS mode (original):
    rnamaps -i rMATS.SE.MATS.JC.txt -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p PTBP1

  VastDB mode (ID lists):
    rnamaps --vastdb_mode \\
      --vastdb_enhanced enhanced.txt --vastdb_silenced silenced.txt \\
      --vastdb_control control.txt --vastdb_constitutive constitutive.txt \\
      --vastdb_annotation EVENT_INFO-hg38.tab \\
      -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p AQR_K562
        """
    )

    # INPUT MODE - mutually exclusive
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-i', '--inputsplice', type=str,
        help='rMATS differential splicing file (rMATS mode)'
    )
    input_group.add_argument(
        '--vastdb_mode', action='store_true',
        help='Use VastDB ID lists mode (requires --vastdb_* arguments)'
    )

    # VASTDB-SPECIFIC ARGUMENTS
    vastdb_group = parser.add_argument_group('VastDB mode options')
    vastdb_group.add_argument(
        '--vastdb_enhanced', help='Enhanced exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_silenced', help='Silenced exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_control', help='Control exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_constitutive', help='Constitutive exon IDs (one per line)')
    vastdb_group.add_argument(
        '--vastdb_annotation', help='VastDB EVENT_INFO file (e.g. EVENT_INFO-hg38.tab)')

    # SHARED REQUIRED ARGUMENTS
    required = parser.add_argument_group('Required arguments (both modes)')
    required.add_argument(
        '-x', '--inputxlsites', type=str, nargs='?',
        help='CLIP crosslinks in BED file format')
    required.add_argument(
        '-f', '--genomefasta', type=str, required=True,
        help='Genome FASTA file (.fa)')
    required.add_argument(
        '-fi', '--fastaindex', type=str, required=True,
        help='Genome FASTA index (.fai)')

    # SHARED OPTIONAL ARGUMENTS
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument(
        '-o', '--outputpath', type=str, default=os.getcwd(), nargs='?',
        help='Output folder [DEFAULT: current directory]')
    optional.add_argument(
        '-w', '--window', type=int, default=300, nargs='?',
        help='Window around splice sites [DEFAULT: 300]')
    optional.add_argument(
        '-s', '--smoothing', type=int, default=15, nargs='?',
        help='Smoothing window [DEFAULT: 15]')
    optional.add_argument(
        '--seed', type=int, default=42,
        help='Random seed for reproducible subsetting [DEFAULT: 42]')
    optional.add_argument(
        '-nc', '--no_constitutive', action='store_true',
        help='Exclude constitutive category')
    optional.add_argument(
        '-ns', '--no_subset', action='store_true',
        help='Disable subsetting of control/constitutive exons')
    optional.add_argument(
        '-ao', '--all_sites', action='store_true',
        help='Include all 6 splice sites (default: 4 core sites)')
    optional.add_argument(
        '-p', '--prefix', type=str,
        help='Prefix for output files')
    optional.add_argument(
        '--hg38_chr_autodetect', action='store_true',
        help='Auto-detect mismatched chromosome naming between the CLIP BED '
             'file (-x) and the exon coordinates, and convert the BED file '
             "using a two-column mapping file (Ensembl <-> GENCODE 'chr' "
             'style). Default mapping file: '
             'test/GRCh38_ensembl2gencode.txt')
    optional.add_argument(
        '--chr_mapping_file', type=str,
        default=os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            'test', 'GRCh38_ensembl2gencode.txt'),
        help='Two-column TSV mapping Ensembl chrom names to GENCODE chrom '
             'names (used with --hg38_chr_autodetect)')

    # PERMUTATION TEST OPTIONS
    perm_group = parser.add_argument_group('Permutation test options')
    perm_group.add_argument(
        '--permute', dest='permute', action='store_true', default=True,
        help='Use label-permutation test for p-values [DEFAULT: on]')
    perm_group.add_argument(
        '--no-permute', dest='permute', action='store_false',
        help="Disable permutation; fall back to per-position Fisher's exact "
             "test (legacy behaviour)")
    perm_group.add_argument(
        '--n_perm', type=int, default=1000,
        help='Number of label permutations [DEFAULT: 1000]')

    # rMATS-SPECIFIC THRESHOLDS
    rmats_group = parser.add_argument_group('rMATS mode thresholds')
    rmats_group.add_argument(
        '-mc', '--minctrl', type=float, default=-0.05, nargs='?',
        help='Minimum dPSI for control events [DEFAULT: -0.05]')
    rmats_group.add_argument(
        '-xc', '--maxctrl', type=float, default=0.05, nargs='?',
        help='Maximum dPSI for control events [DEFAULT: 0.05]')
    rmats_group.add_argument(
        '-xi', '--maxincl', type=float, default=0.9, nargs='?',
        help='Maximum PSI for control (above = constitutive) [DEFAULT: 0.9]')
    rmats_group.add_argument(
        '-xf', '--maxfdr', type=float, default=0.1, nargs='?',
        help='Maximum FDR for regulated events [DEFAULT: 0.1]')
    rmats_group.add_argument(
        '-xe', '--maxenh', type=float, default=-0.05, nargs='?',
        help='Maximum dPSI for enhanced exons [DEFAULT: -0.05]')
    rmats_group.add_argument(
        '-ms', '--minsil', type=float, default=0.05, nargs='?',
        help='Minimum dPSI for silenced exons [DEFAULT: 0.05]')

    # MULTIVALENCY (rMATS mode)
    mv_group = parser.add_argument_group('Multivalency analysis')
    mv_group.add_argument(
        '-v', '--multivalency', action='store_true',
        help='Run multivalency analysis (requires germs.R)')
    mv_group.add_argument(
        '-g', '--germsdir', type=str, default=os.getcwd(), nargs='?',
        help='Directory containing germs.R [DEFAULT: current directory]')

    args = parser.parse_args()

    # Validate VastDB mode requirements
    if args.vastdb_mode:
        if not args.vastdb_annotation:
            parser.error("--vastdb_mode requires --vastdb_annotation "
                         "(EVENT_INFO file)")
        if not any([args.vastdb_enhanced, args.vastdb_silenced,
                    args.vastdb_control, args.vastdb_constitutive]):
            parser.error("--vastdb_mode requires at least one "
                         "--vastdb_* ID list file")

    return args

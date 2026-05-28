"""Argparse-based command line interface."""

import argparse
import os


def cli():
    parser = argparse.ArgumentParser(
        prog='rnamaps',
        description='Plot CLIP crosslinks around regulated exons to study '
                    'position-dependent impact on pre-mRNA splicing.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Input modes:

        rMATS mode (default):
            rnamaps -i rMATS.SE.MATS.JC.txt -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p PTBP1

        VastDB mode (ID lists):
            rnamaps --vastdb_mode \
            --vastdb_enhanced enhanced.txt --vastdb_silenced silenced.txt \
            --vastdb_control control.txt --vastdb_constitutive constitutive.txt \
            --vastdb_annotation EVENT_INFO-hg38.tab \
            -x CLIP.bed -f hg38.fa -fi hg38.fa.fai -o output -p AQR_K562
        """
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-i', '--inputsplice', type=str,
        help='rMATS differential splicing file (rMATS mode)'
    )
    input_group.add_argument(
        '--vastdb_mode', action='store_true',
        help='Use VastDB ID lists mode (requires --vastdb_* arguments)'
    )
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument(
        '--y_axis', type=str, default='log10p', choices=['log10p', 'zscore'],
        help="Y-axis for legacy permutation_z RNA map plots: 'log10p' "
             "(default) for signed -log10(p), 'zscore' for signed "
             "permutation z-score. Ignored for other --enrichment methods.")
    optional.add_argument(
        '--enrichment', type=str, nargs='+', default=None,
        choices=['permutation_z', 'fisher', 'bootstrap_contrast',
                 'cluster_perm', 'roc_auc'],
        help="One or more enrichment methods to run. "
             "[DEFAULT: bootstrap_contrast]. Honours legacy --no-permute "
             "as 'fisher' when --enrichment is not given.")
    optional.add_argument(
        '--binarise', dest='binarise', action='store_true', default=True,
        help="Threshold the per-exon coverage matrix to 0/1 before "
             "running enrichment (CLIP-style 'did this exon have any "
             "signal at this base?'). [DEFAULT: on]")
    optional.add_argument(
        '--no-binarise', dest='binarise', action='store_false',
        help="Keep raw continuous values in the per-exon coverage "
             "matrix. Use this when the -x input is a continuous "
             "signal (e.g. AI prediction scores or a density track) "
             "and the per-(exon, position) magnitude is itself "
             "informative. The test statistics become differences/"
             "ratios of mean signal levels rather than positive-exon "
             "fractions. Note: with --enrichment bootstrap_contrast "
             "--shrinkage magnitude, the default --shrinkage_scale "
             "(0.05) assumes a [0, 1] rate axis -- override it to "
             "match the scale of your continuous signal.")

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

    # CONTROL-SET HYGIENE OPTIONS
    ctrl_group = parser.add_argument_group('Control-set hygiene options')
    ctrl_group.add_argument(
        '--control_set', type=str, default='default',
        choices=['default', 'strict', 'constitutive_only'],
        help="Control-set mode. 'default' uses the current control "
             "criteria. 'strict' tightens to |dPSI|<--control_max_dpsi "
             "AND FDR>--control_min_fdr. 'constitutive_only' drops the "
             "control category and relabels constitutive -> control."
    )
    ctrl_group.add_argument(
        '--control_max_dpsi', type=float, default=0.01,
        help="Strict mode: max |dPSI| for control exons [DEFAULT: 0.01]"
    )
    ctrl_group.add_argument(
        '--control_min_fdr', type=float, default=0.5,
        help="Strict mode: min FDR for control exons [DEFAULT: 0.5]"
    )

    # EXPRESSION-MATCHING OPTIONS
    expr_group = parser.add_argument_group('Expression-matching options')
    expr_group.add_argument(
        '--gene_tpm', type=str, default=None,
        help='Optional 2-column gene TPM table (gene_id,tpm). '
             'When provided, control/constitutive sets are expression-matched '
             'to regulated exons.'
    )
    expr_group.add_argument(
        '--tpm_n_bins', type=int, default=10,
        help='Number of quantile bins for regulated log-TPM distribution '
             '[DEFAULT: 10]'
    )
    expr_group.add_argument(
        '--tpm_pseudocount', type=float, default=1.0,
        help='Pseudocount for log10(TPM + pseudocount) transform '
             '[DEFAULT: 1.0]'
    )
    expr_group.add_argument(
        '--tpm_min_tpm', type=float, default=0.0,
        help='Minimum TPM required before matching [DEFAULT: 0.0]'
    )
    expr_group.add_argument(
        '--no_match_constitutive', action='store_true',
        help='By default, expression matching applies to control and '
             'constitutive sets. Use this flag to leave constitutive '
             'untouched.'
    )

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
        help='Number of label permutations [DEFAULT: 1000] (used by '
             'permutation_z and cluster_perm)')

    # BOOTSTRAP CONTRAST OPTIONS
    boot_group = parser.add_argument_group(
        'Bootstrap contrast options (--enrichment bootstrap_contrast)')
    boot_group.add_argument(
        '--n_boot', type=int, default=1000,
        help='Bootstrap iterations [DEFAULT: 1000]')
    boot_group.add_argument(
        '--bootstrap_control_fixed', action='store_true',
        help='Treat control mean as a constant (skip resampling control). '
             'Equivalent up to negligible variance when n_ctrl >> n_c '
             'and 5-10x faster.')
    boot_group.add_argument(
        '--shrinkage', type=str, default='none',
        choices=['none', 'magnitude', 'pseudocount'],
        help="How log2fc / log_odds_ratio are regularised. 'none' "
             "(default): no shrinkage. log2fc may produce +-inf "
             "where one rate is exactly zero; log_odds_ratio clips "
             "rates into [eps_safe, 1-eps_safe] only to keep the "
             "logit finite. 'magnitude': both log2fc and "
             "log_odds_ratio are shrunk toward 0 by a weight that "
             "depends on the LARGER of the two rates -- not on the "
             "number of exons. When both rates are small in absolute "
             "terms (max(cat, ctrl) below the shrinkage scale), the "
             "contrasts are pulled toward zero; when at least one "
             "rate is well above the scale, the raw contrast is "
             "preserved. 'pseudocount': legacy additive epsilon. "
             "Note: the previous default 'magnitude' assumed a "
             "[0,1] rate axis (binarised CLIP); for continuous "
             "inputs via --no-binarise it produced confusing "
             "shrinkage, hence the new 'none' default.")
    boot_group.add_argument(
        '--shrinkage_scale', type=float, default=None,
        help='Only used with --shrinkage magnitude. Sets the rate '
             'scale tau at which shrinkage transitions from heavy to '
             'light. Default (None): tau = 0.05 (5%% rate). Below this '
             'rate, log2fc / log_odds_ratio are pulled toward zero; '
             'well above it, the raw contrast is preserved. Pass a '
             'smaller value (e.g. 0.01) for gentler shrinkage that '
             'only bites at very small rates, or a larger value for '
             'more aggressive shrinkage even at moderate rates.')
    boot_group.add_argument(
        '--pseudocount', type=float, default=None,
        help='Only used with --shrinkage pseudocount. Override adaptive '
             'log2FC pseudocount with a fixed value. When unset, '
             'eps = max(1e-3, --pseudocount_frac * '
             'median(mean_cov_ctrl over region)).')
    boot_group.add_argument(
        '--pseudocount_frac', type=float, default=0.01,
        help='Only used with --shrinkage pseudocount. Adaptive log2FC '
             'pseudocount fraction [DEFAULT: 0.01]')

    # CLUSTER-PERMUTATION OPTIONS
    cl_group = parser.add_argument_group(
        'Cluster-permutation options (--enrichment cluster_perm)')
    cl_group.add_argument(
        '--cluster_thresh', type=float, default=2.0,
        help='Cluster-defining |t| threshold [DEFAULT: 2.0]')

    # BED-SCORE / ROC-AUC OPTIONS
    score_group = parser.add_argument_group(
        'BED score handling and ROC/AUC options (--enrichment roc_auc)')
    score_group.add_argument(
        '--xl_score', type=str, nargs='+', default=['ignore'],
        choices=['ignore', 'raw', 'per_transcript_zscore'],
        help="How to use BED column 5 of -x. Accepts one or more "
             "modes; passing several runs each requested enrichment "
             "method once per mode and writes per-mode output files. "
             "The rMATS categorisation, exon-length plot, heatmap, "
             "and totalExonsCovered table are xl_score-invariant and "
             "are written only once. "
             "'ignore' (default) counts each overlap as 1, matching "
             "the legacy behaviour. 'raw' uses the BED score as-is "
             "(good for eCLIP integer counts or any pre-normalised "
             "track). 'per_transcript_zscore' z-scores each "
             "transcript's score row before aggregation (good for AI "
             "prediction tracks where absolute score magnitudes are "
             "not comparable across transcripts). The per-transcript "
             "mode uses BED column 4 (name) as the transcript "
             "identifier; rows with placeholder '.' names are kept "
             "but counted as their own singleton transcripts.")
    score_group.add_argument(
        '--roc_aggregator', type=str, default='both',
        choices=['mean', 'max', 'both'],
        help="How to aggregate per-exon signal over the window for the "
             "per-region ROC curve. 'mean' (or equivalently 'sum' for a "
             "fixed-width window) rewards broad coverage; 'max' rewards "
             "sharp localised peaks. 'both' computes mean and max, "
             "writes both to TSV, and plots mean. [DEFAULT: both]")
    score_group.add_argument(
        '--roc_n_perm', type=int, default=0,
        help="Number of label permutations for AUC null distribution. "
             "0 (default) skips the null and reports only the observed "
             "AUC; >0 also reports a permutation p-value per position "
             "and per region.")

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

    # Resolve --enrichment with backward-compat for legacy --permute /
    # --no-permute. Explicit --enrichment wins. Otherwise, --no-permute
    # maps to 'fisher'; default is 'bootstrap_contrast'.
    if args.enrichment is None:
        if not args.permute:
            args.enrichment = ['fisher']
        else:
            args.enrichment = ['bootstrap_contrast']

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

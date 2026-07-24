"""Top-level RNA map pipeline that orchestrates loading, processing, and plotting."""

import logging
import os
import sys
import traceback

import numpy as np
import pandas as pd
import pybedtools as pbt

from rnamaps.coverage import detect_nontrivial_bed_scores, get_coverage_plot
from rnamaps.enrichment import EnrichmentResult
from rnamaps.enrichment import bin_tpr as enrich_bin_tpr
from rnamaps.enrichment import bootstrap_contrast as enrich_bootstrap
from rnamaps.enrichment import cluster_perm as enrich_cluster
from rnamaps.enrichment import roc_auc as enrich_roc_auc
from rnamaps.expression_matching import (
    attach_tpm_to_exons,
    load_tpm_table,
    match_controls_by_expression,
)
from rnamaps.io_rmats import load_rmats_data
from rnamaps.io_vastdb import load_vastdb_data
from rnamaps.logging_utils import log_runtime, setup_logging
from rnamaps.multivalency import plot_multivalency
from rnamaps.permutation import compute_permutation_pvalues
from rnamaps.plots import (
    plot_exon_lengths,
    plot_heatmap,
    plot_rna_map,
    plot_roc_curves,
)
from rnamaps.preprocessing import (
    apply_control_set,
    apply_subsetting,
    autodetect_and_convert_bed_chroms,
    autodetect_and_convert_df_chroms,
    decompress_if_gz,
    get_ss_bed,
)


def _run_method(method, region_cov, region_fisher_linegraph,
                exon_categories, region_label, args, rng, smoothing):
    """Dispatch a single enrichment method on one splice-site region.

    Returns
    -------
    EnrichmentResult or None
    """
    binarise = getattr(args, 'binarise', True)
    if method == 'bootstrap_contrast':
        return enrich_bootstrap.compute(
            region_cov, exon_categories, region_label,
            rng=rng,
            n_boot=args.n_boot,
            shrinkage=getattr(args, 'shrinkage', 'magnitude'),
            shrinkage_scale=getattr(args, 'shrinkage_scale', None),
            pseudocount=args.pseudocount,
            pseudocount_frac=args.pseudocount_frac,
            bootstrap_control_fixed=args.bootstrap_control_fixed,
            smoothing=smoothing,
            binarise=binarise,
        )
    if method == 'cluster_perm':
        return enrich_cluster.compute(
            region_cov, exon_categories, region_label,
            rng=rng,
            n_perm=args.n_perm,
            cluster_thresh=args.cluster_thresh,
            binarise=binarise,
        )
    if method == 'roc_auc':
        return enrich_roc_auc.compute(
            region_cov, exon_categories, region_label,
            rng=rng,
            roc_aggregator=getattr(args, 'roc_aggregator', 'both'),
            n_perm=getattr(args, 'roc_n_perm', 0),
            smoothing=smoothing,
            binarise=False,
        )
    if method == 'bin_tpr':
        return enrich_bin_tpr.compute(
            region_cov, exon_categories, region_label,
            bin_size=getattr(args, 'bin_tpr_size', 50),
        )
    if method == 'permutation_z':
        plot_df, clusters_df = compute_permutation_pvalues(
            region_cov, exon_categories, region_label,
            n_perm=args.n_perm, smoothing=smoothing, rng=rng,
            binarise=binarise,
        )
        y_col = ('zscore' if getattr(args, 'y_axis', 'log10p') == 'zscore'
                 else '-log10pvalue_smoothed')
        ylabel = ('signed permutation z-score vs control'
                  if y_col == 'zscore'
                  else 'signed -log10(empirical p) vs control')
        return EnrichmentResult(
            plot_df=plot_df, clusters_df=clusters_df,
            plot_kind='line', y_columns=(y_col,),
            method_name='permutation_z', ylabel=ylabel,
        )
    if method == 'fisher':
        df = region_fisher_linegraph
        if df is None or df.empty:
            return None
        return EnrichmentResult(
            plot_df=df, clusters_df=pd.DataFrame(),
            plot_kind='line', y_columns=('-log10pvalue_smoothed',),
            method_name='fisher',
            ylabel='-log10(p value) enrichment / control',
        )
    raise ValueError(f"Unknown enrichment method: {method!r}")


def _plot_method(method, plot_df, clusters_df, exon_categories,
                 original_counts, window, args, output_dir, FILEname,
                 roc_curves_df=None, region_auc_df=None,
                 xl_score_mode=None):
    """Render one PDF (or two for bootstrap_contrast) for a method's
    accumulated per-region results.

    ``xl_score_mode`` is the BED-score handling used to build the
    underlying coverage matrices; it's woven into the output PDF
    filenames only when the pipeline is sweeping multiple modes (so
    single-mode runs keep their backward-compatible names)."""
    suffix = f"_xlscore-{xl_score_mode}" if xl_score_mode else ""
    if method == 'bootstrap_contrast':
        smooth_str = (f", smooth={args.smoothing}"
                      if getattr(args, 'smoothing', 1)
                      and args.smoothing > 1 else "")
        ctrl_str = (' [ctrl fixed]'
                    if args.bootstrap_control_fixed else '')
        shrinkage_mode = getattr(args, 'shrinkage', 'magnitude')
        shrinkage_str = f", shrinkage={shrinkage_mode}"
        subtitle = (f"Bootstrap contrast (B={args.n_boot}"
                    f"{smooth_str}{shrinkage_str}){ctrl_str}")
        for y_col, ylab in [
            ('delta',
             'fraction of category exons positive - control fraction'),
            ('log2fc',
             'log2 fold change vs control'),
            ('log_odds_ratio',
             'log odds ratio (logit(cat) - logit(ctrl))'),
        ]:
            plot_rna_map(
                plot_df, exon_categories, original_counts,
                window, args.all_sites, output_dir, FILEname,
                plot_kind='ribbon', y_col=y_col,
                ci_cols=(f'{y_col}_lo', f'{y_col}_hi'),
                method_name=f"{y_col}{suffix}",
                ylabel=ylab, subtitle=subtitle,
            )
        return
    if method == 'cluster_perm':
        subtitle = (f"Cluster-mass permutation (B={args.n_perm}, "
                    f"|t|>{args.cluster_thresh}); bars at p<=0.05")
        plot_rna_map(
            plot_df, exon_categories, original_counts,
            window, args.all_sites, output_dir, FILEname,
            plot_kind='clusters', y_col='t_obs',
            clusters_df=clusters_df,
            method_name=f'cluster_perm{suffix}',
            ylabel='Welch t (cluster permutation)',
            subtitle=subtitle,
        )
        return
    if method == 'permutation_z':
        plot_rna_map(
            plot_df, exon_categories, original_counts,
            window, args.all_sites, output_dir, FILEname,
            pvalue_method='permutation', n_perm=args.n_perm,
            y_axis=getattr(args, 'y_axis', 'log10p'),
            plot_kind='line', method_name=f'permutation_z{suffix}',
        )
        return
    if method == 'fisher':
        plot_rna_map(
            plot_df, exon_categories, original_counts,
            window, args.all_sites, output_dir, FILEname,
            pvalue_method='fisher', plot_kind='line',
            method_name=f'fisher{suffix}',
        )
        return
    if method == 'bin_tpr':
        # Tabular-only output; the TSV is written by the caller.
        return
    if method == 'roc_auc':
        score_mode = xl_score_mode or 'ignore'
        n_perm = getattr(args, 'roc_n_perm', 0)
        agg_choice = getattr(args, 'roc_aggregator', 'both')
        subtitle = (
            f"Per-position AUC (xl_score={score_mode}, "
            f"aggregator={agg_choice}"
            + (f", n_perm={n_perm}" if n_perm > 0 else "")
            + ")"
        )
        # 1. Per-position AUC line plot.
        plot_rna_map(
            plot_df, exon_categories, original_counts,
            window, args.all_sites, output_dir, FILEname,
            plot_kind='line', y_col='auc_signed_smoothed',
            method_name=f'roc_auc{suffix}',
            ylabel='signed AUC (2*(AUC - 0.5)) vs control',
            subtitle=subtitle,
        )
        # 2. Per-region ROC curves.
        if roc_curves_df is not None and not roc_curves_df.empty:
            plot_roc_curves(
                roc_curves_df, region_auc_df, args.all_sites,
                output_dir, FILEname,
                aggregators_to_plot=(['mean', 'max']
                                     if agg_choice == 'both'
                                     else [agg_choice]),
                xl_score_mode=xl_score_mode,
            )
        return
    raise ValueError(f"Unknown enrichment method: {method!r}")


def run_rna_map(args):
    """
    Main RNA map pipeline. Handles both input modes with shared downstream logic.
    """
    output_dir = args.outputpath
    os.makedirs(output_dir, exist_ok=True)

    log_filename, start_time, logger = setup_logging(output_dir)
    logging.info(f"Log file created: {log_filename}")
    logging.info(f"Arguments: {args}")

    try:
        # Set random seed for reproducibility (legacy global seed for any
        # downstream code that still uses np.random.*) and a Generator for
        # the permutation test.
        np.random.seed(args.seed)
        rng = np.random.default_rng(args.seed)
        logging.info(f"Random seed set to {args.seed}")

        # Load chromosome list
        df_fai = pd.read_csv(args.fastaindex, sep='\t', header=None)
        chroms = set(df_fai[0].values)

        # ==============================================================
        # MODE SELECTION: Two tracks, one output format
        # ==============================================================

        if args.inputsplice:
            # ----- TRACK 1: rMATS -----
            input_mode = 'rmats'
            df_rmats = load_rmats_data(
                args.inputsplice,
                args.minctrl, args.maxctrl, args.maxincl,
                args.maxfdr, args.maxenh, args.minsil,
                chroms, args.no_constitutive
            )

            if args.prefix:
                FILEname = (args.prefix + "_" +
                            args.inputsplice.split('/')[-1]
                            .replace('.txt', '').replace('.gz', ''))
            else:
                FILEname = (args.inputsplice.split('/')[-1]
                            .replace('.txt', '').replace('.gz', ''))

        else:
            # ----- TRACK 2: VastDB -----
            input_mode = 'vastdb'
            df_rmats = load_vastdb_data(
                args.vastdb_enhanced,
                args.vastdb_silenced,
                args.vastdb_control,
                args.vastdb_constitutive,
                args.vastdb_annotation,
                chroms
            )

            if args.prefix:
                FILEname = args.prefix
            else:
                FILEname = "VastDB_RNAmap"

        # ==============================================================
        # SHARED PIPELINE: Same for both modes from here on
        # ==============================================================

        # Auto-convert exon chrom naming to match fai if requested
        if getattr(args, 'hg38_chr_autodetect', False):
            df_rmats = autodetect_and_convert_df_chroms(
                df_rmats, chroms, args.chr_mapping_file,
                chr_col='chr', label='exon')

        # Filter to valid chromosomes
        df_rmats = df_rmats[df_rmats['chr'].isin(chroms)]

        # Remove exons with missing flanking coordinates
        logging.info("\nFiltering exons with complete flanking coordinates...")
        before = len(df_rmats)
        df_rmats = df_rmats.dropna(
            subset=['upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
        )
        after = len(df_rmats)
        logging.info(f"Removed {before - after} exons with missing flanking coordinates")
        logging.info(f"Remaining: {after} exons")

        # Apply control-set hygiene mode (no-op for 'default').
        df_rmats = apply_control_set(
            df_rmats,
            mode=getattr(args, 'control_set', 'default'),
            control_max_dpsi=getattr(args, 'control_max_dpsi', 0.01),
            control_min_fdr=getattr(args, 'control_min_fdr', 0.5),
        )

        if getattr(args, 'gene_tpm', None):
            logging.info("\nApplying expression-aware control matching...")
            tpm_df = load_tpm_table(args.gene_tpm)
            df_rmats, attach_info = attach_tpm_to_exons(df_rmats, tpm_df)
            logging.info(f"TPM attachment summary: {attach_info}")
            df_rmats, match_summary = match_controls_by_expression(
                df_rmats,
                n_bins=getattr(args, 'tpm_n_bins', 10),
                pseudocount=getattr(args, 'tpm_pseudocount', 1.0),
                min_tpm=getattr(args, 'tpm_min_tpm', 0.0),
                also_constitutive=not getattr(
                    args, 'no_match_constitutive', False
                ),
                rng=rng,
            )
            logging.info(f"Expression-matching summary: {match_summary}")

        exon_categories = df_rmats.groupby('category').size()
        logging.info("\nExons in each category:")
        logging.info(exon_categories)

        # Validate categories
        if "control" not in exon_categories or exon_categories.loc["control"] == 0:
            logging.error(
                "No control exons found! If using --gene_tpm, relax "
                "--tpm_min_tpm and/or reduce --tpm_n_bins."
            )
            sys.exit(1)

        if ("enhanced" not in exon_categories
                and "silenced" not in exon_categories):
            logging.error("No regulated exons found!")
            sys.exit(1)

        # Apply subsetting. Only the legacy 'fisher' method benefits
        # from subsetting; bootstrap_contrast / cluster_perm /
        # permutation_z all handle unequal n correctly. Auto-disable
        # subsetting whenever any non-fisher method is selected so the
        # legend doesn't show stale "subset from N" annotations.
        enrichment_methods = list(getattr(args, 'enrichment',
                                          ['bootstrap_contrast']))
        needs_full = any(m != 'fisher' for m in enrichment_methods)
        if not args.no_subset and not needs_full:
            df_rmats, original_counts = apply_subsetting(
                df_rmats, args.no_constitutive
            )
        else:
            if needs_full and not args.no_subset:
                logging.info(
                    "Subsetting auto-disabled because --enrichment "
                    "includes a non-fisher method "
                    f"({', '.join(enrichment_methods)})."
                )
            else:
                logging.info("Subsetting disabled (--no_subset flag)")
            category_counts = df_rmats['category'].value_counts()
            original_counts = {cat: count for cat, count in category_counts.items()}

        exon_categories = df_rmats.groupby('category').size()

        # Save categorised exons
        if input_mode == 'rmats':
            save_cols = ['chr', 'exonStart_0base', 'exonEnd', 'strand', 'category',
                         'FDR', 'dPSI', 'maxPSI', 'GeneID', 'geneSymbol',
                         'gene_id', 'tpm',
                         'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
            save_cols = [c for c in save_cols if c in df_rmats.columns]
            suffix = '_RMATS_with_categories.tsv'
        else:
            save_cols = ['chr', 'exonStart_0base', 'exonEnd', 'strand', 'category',
                         'EVENT', 'GENE', 'gene_id', 'tpm',
                         'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
            save_cols = [c for c in save_cols if c in df_rmats.columns]
            suffix = '_VastDB_with_categories.tsv'

        df_rmats[save_cols].to_csv(
            f'{output_dir}/{FILEname}{suffix}', sep="\t", index=False
        )
        logging.info(f"Saved categorised exons to {output_dir}/{FILEname}{suffix}")

        # Exon length plots
        plot_exon_lengths(df_rmats.copy(), output_dir, FILEname)

        # ==============================================================
        # CREATE BED FILES FOR 6 SPLICE SITE REGIONS
        # ==============================================================
        # After strand correction at load time, upstream/downstream are
        # always in transcript order. Same-exon pairing is correct for
        # both modes.
        # ==============================================================
        logging.info("\n" + "=" * 60)
        logging.info("CREATING BED FILES FOR SPLICE SITES")
        logging.info("=" * 60)

        middle_3ss_bed = get_ss_bed(df_rmats, 'exonStart_0base', 'exonEnd')
        middle_5ss_bed = get_ss_bed(df_rmats, 'exonEnd', 'exonStart_0base')
        downstream_3ss_bed = get_ss_bed(df_rmats, 'downstreamES', 'downstreamEE')
        upstream_5ss_bed = get_ss_bed(df_rmats, 'upstreamEE', 'upstreamES')

        downstream_5ss_bed = None
        upstream_3ss_bed = None
        if args.all_sites:
            downstream_5ss_bed = get_ss_bed(df_rmats, 'downstreamEE', 'downstreamES')
            upstream_3ss_bed = get_ss_bed(df_rmats, 'upstreamES', 'upstreamEE')

        # ==============================================================
        # CALCULATE COVERAGE
        # ==============================================================
        if args.inputxlsites is not None:
            logging.info("\n" + "=" * 60)
            logging.info("CALCULATING COVERAGE")
            logging.info("=" * 60)

            fai = args.fastaindex
            xl_bed = args.inputxlsites
            window = args.window
            smoothing = args.smoothing

            xl_bed = decompress_if_gz(xl_bed, output_dir)

            if getattr(args, 'hg38_chr_autodetect', False):
                xl_bed = autodetect_and_convert_bed_chroms(
                    xl_bed, chroms, args.chr_mapping_file, output_dir)

            # Normalise --xl_score (either a single string -- legacy
            # access from tests/notebooks -- or a list-of-strings from
            # the CLI). Single 'ignore' mode keeps backward-compatible
            # filenames; multi-mode runs add a per-mode suffix.
            xl_score_arg = getattr(args, 'xl_score', ['ignore'])
            if isinstance(xl_score_arg, str):
                xl_score_modes = [xl_score_arg]
            else:
                xl_score_modes = list(xl_score_arg) or ['ignore']
            multi_mode = len(xl_score_modes) > 1

            if 'ignore' in xl_score_modes and detect_nontrivial_bed_scores(
                xl_bed
            ):
                logging.warning(
                    "BED column 5 of %s contains non-trivial scores. "
                    "--xl_score includes 'ignore', which treats them as "
                    "presence indicators only. Add (or switch to) --xl_score "
                    "raw / per_transcript_zscore "
                    "to use the score values.",
                    xl_bed,
                )
            logging.info(f"BED score modes: {xl_score_modes}")

            # ----------------------------------------------------------
            # 1. Heatmap and totalExonsCovered are xl_score-invariant.
            #    Compute them once from a binarised "presence" matrix
            #    (xl_score=ignore). Skip if user explicitly only asked
            #    for score-aware modes; the heatmap is a "which exons
            #    have any signal" view, so ignore-mode is always the
            #    right basis.
            # ----------------------------------------------------------
            logging.info("Computing presence-only coverage for heatmap...")
            heatmap_region_covs = []
            for region_label_, region_bed_ in [
                ('upstream_5ss', upstream_5ss_bed),
                ('middle_3ss', middle_3ss_bed),
                ('middle_5ss', middle_5ss_bed),
                ('downstream_3ss', downstream_3ss_bed),
            ] + ([
                ('downstream_5ss', downstream_5ss_bed),
                ('upstream_3ss', upstream_3ss_bed),
            ] if args.all_sites else []):
                _, rc_ = get_coverage_plot(
                    xl_bed, region_bed_, fai, window, exon_categories,
                    region_label_, smoothing, score_mode='ignore',
                )
                heatmap_region_covs.append(rc_)
            if getattr(args, 'dump_region_matrix', False):
                for rc_ in heatmap_region_covs:
                    npz_path = (
                        f"{output_dir}/{FILEname}_{rc_.label}"
                        f"_region_matrix.npz"
                    )
                    rc_.save_npz(npz_path)
                    logging.info(
                        f"Dumped region matrix {rc_.label} "
                        f"({rc_.n_exons}×{rc_.n_positions}) to {npz_path}"
                    )
            plot_heatmap(heatmap_region_covs, exon_categories, window,
                         args.all_sites, output_dir, FILEname)

            # ----------------------------------------------------------
            # 2. Per xl_score mode: compute coverage matrices, run all
            #    requested enrichment methods, emit per-method outputs
            #    (with the score-mode suffix when sweeping multiple).
            # ----------------------------------------------------------
            for score_mode in xl_score_modes:
                file_suffix = (f"_xlscore-{score_mode}" if multi_mode
                               else "")
                logging.info("\n" + "=" * 60)
                logging.info(
                    f"COVERAGE / ENRICHMENT (--xl_score {score_mode})"
                )
                logging.info("=" * 60)

                try:
                    if score_mode == 'ignore':
                        # Reuse the presence-only matrices already built
                        # for the heatmap so we don't pay bedtools twice.
                        region_covs_mode = list(heatmap_region_covs)
                        fisher_linegraphs_mode = []
                        for rc_ in region_covs_mode:
                            # Re-run the fast Fisher-linegraph builder
                            # from the matrix (no bedtools needed).
                            from rnamaps.coverage import (
                                _fisher_linegraph_from_matrix,
                            )
                            fisher_linegraphs_mode.append(
                                _fisher_linegraph_from_matrix(rc_, smoothing)
                            )
                    else:
                        region_covs_mode = []
                        fisher_linegraphs_mode = []
                        region_specs = [
                            ('upstream_5ss', upstream_5ss_bed),
                            ('middle_3ss', middle_3ss_bed),
                            ('middle_5ss', middle_5ss_bed),
                            ('downstream_3ss', downstream_3ss_bed),
                        ]
                        if args.all_sites:
                            region_specs.extend([
                                ('downstream_5ss', downstream_5ss_bed),
                                ('upstream_3ss', upstream_3ss_bed),
                            ])
                        for region_label_, region_bed_ in region_specs:
                            lg_, rc_ = get_coverage_plot(
                                xl_bed, region_bed_, fai, window,
                                exon_categories, region_label_, smoothing,
                                score_mode=score_mode,
                            )
                            fisher_linegraphs_mode.append(lg_)
                            region_covs_mode.append(rc_)
                        if getattr(args, 'dump_region_matrix', False):
                            for rc_ in region_covs_mode:
                                npz_path = (
                                    f"{output_dir}/{FILEname}_{rc_.label}"
                                    f"_xlscore-{score_mode}_region_matrix.npz"
                                )
                                rc_.save_npz(npz_path)
                                logging.info(
                                    f"Dumped region matrix "
                                    f"{rc_.label}/{score_mode} "
                                    f"({rc_.n_exons}×{rc_.n_positions}) "
                                    f"to {npz_path}"
                                )

                    # method -> list of EnrichmentResult per region.
                    results_by_method = {m: [] for m in enrichment_methods}
                    for region_cov, fisher_lg in zip(
                        region_covs_mode, fisher_linegraphs_mode
                    ):
                        if region_cov.n_exons == 0:
                            continue
                        region_label = region_cov.label
                        for method in enrichment_methods:
                            res = _run_method(
                                method, region_cov, fisher_lg,
                                exon_categories, region_label, args, rng,
                                smoothing,
                            )
                            if res is not None:
                                results_by_method[method].append(res)

                    for method in enrichment_methods:
                        method_results = results_by_method[method]
                        if not method_results:
                            logging.warning(
                                f"[{method}/{score_mode}] No results to plot."
                            )
                            continue

                        method_plot_df = pd.concat(
                            [r.plot_df for r in method_results],
                            ignore_index=True,
                        )
                        cluster_frames = [
                            r.clusters_df for r in method_results
                            if r.clusters_df is not None
                            and not r.clusters_df.empty
                        ]
                        method_clusters_df = (
                            pd.concat(cluster_frames, ignore_index=True)
                            if cluster_frames else pd.DataFrame()
                        )

                        roc_curves_df = None
                        region_auc_df = None
                        if any(getattr(r, 'extras', None)
                               for r in method_results):
                            roc_frames = [
                                r.extras.get('roc_curves_df')
                                for r in method_results
                                if r.extras
                                and r.extras.get('roc_curves_df') is not None
                                and not r.extras['roc_curves_df'].empty
                            ]
                            if roc_frames:
                                roc_curves_df = pd.concat(
                                    roc_frames, ignore_index=True
                                )
                            auc_frames = [
                                r.extras.get('region_auc_df')
                                for r in method_results
                                if r.extras
                                and r.extras.get('region_auc_df') is not None
                                and not r.extras['region_auc_df'].empty
                            ]
                            if auc_frames:
                                region_auc_df = pd.concat(
                                    auc_frames, ignore_index=True
                                )

                        method_plot_df.to_csv(
                            f'{output_dir}/{FILEname}_RNAmap_{method}'
                            f'{file_suffix}.tsv',
                            sep="\t", index=False,
                        )
                        if not method_clusters_df.empty:
                            method_clusters_df.to_csv(
                                f'{output_dir}/{FILEname}_RNAmap_{method}'
                                f'{file_suffix}_clusters.tsv',
                                sep="\t", index=False,
                            )
                        if roc_curves_df is not None and not roc_curves_df.empty:
                            roc_curves_df.to_csv(
                                f'{output_dir}/{FILEname}_RNAmap_{method}'
                                f'{file_suffix}_roc_curves.tsv',
                                sep="\t", index=False,
                            )
                        if region_auc_df is not None and not region_auc_df.empty:
                            region_auc_df.to_csv(
                                f'{output_dir}/{FILEname}_RNAmap_{method}'
                                f'{file_suffix}_region_auc.tsv',
                                sep="\t", index=False,
                            )

                        _plot_method(
                            method, method_plot_df, method_clusters_df,
                            exon_categories, original_counts,
                            window, args, output_dir, FILEname,
                            roc_curves_df=roc_curves_df,
                            region_auc_df=region_auc_df,
                            xl_score_mode=(score_mode if multi_mode else None),
                        )
                except Exception:
                    # Log full traceback to the .log file so future
                    # silent crashes show up here instead of only on
                    # stderr. Skip to the next mode rather than killing
                    # the whole sweep.
                    logging.error(
                        "Unhandled exception while processing "
                        f"--xl_score {score_mode}; skipping this mode. "
                        "Traceback:\n%s",
                        traceback.format_exc(),
                    )
                finally:
                    # Release pybedtools temp files between modes so
                    # /tmp doesn't fill up across a multi-mode sweep
                    # (each bedtools intersect can write a multi-GB
                    # temp file that otherwise lives until process
                    # exit).
                    try:
                        pbt.cleanup(remove_all=False)
                    except Exception:
                        logging.warning(
                            "pybedtools.cleanup() raised; continuing.",
                            exc_info=True,
                        )

        # ==============================================================
        # MULTIVALENCY (optional, requires germs.R)
        # ==============================================================
        if hasattr(args, 'multivalency') and args.multivalency:
            logging.info("\n" + "=" * 60)
            logging.info("MULTIVALENCY ANALYSIS")
            logging.info("=" * 60)

            plot_multivalency(
                middle_3ss_bed, middle_5ss_bed,
                downstream_3ss_bed, upstream_5ss_bed,
                downstream_5ss_bed, upstream_3ss_bed,
                args.fastaindex, args.window, args.genomefasta,
                output_dir, FILEname, args.germsdir,
                args.all_sites, exon_categories, original_counts
            )

        logging.info("\n" + "=" * 60)
        logging.info("SCRIPT COMPLETED SUCCESSFULLY")
        logging.info("=" * 60)

    except Exception:
        # Mirror the traceback into the .log file before re-raising so
        # that crashes are not silent in the per-run log (previously
        # the .log only captured logging.* calls and the traceback
        # went to stderr, leaving no record on disk if stderr was lost).
        logging.error(
            "Unhandled exception in run_rna_map; aborting. "
            "Traceback:\n%s",
            traceback.format_exc(),
        )
        raise
    finally:
        log_runtime(start_time, logger)
        for handler in logger.handlers:
            handler.flush()
        logging.shutdown()

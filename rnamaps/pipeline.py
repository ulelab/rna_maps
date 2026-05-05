"""Top-level RNA map pipeline that orchestrates loading, processing, and plotting."""

import logging
import os
import sys

import numpy as np
import pandas as pd

from rnamaps.coverage import get_coverage_plot
from rnamaps.enrichment import EnrichmentResult
from rnamaps.enrichment import bootstrap_contrast as enrich_bootstrap
from rnamaps.enrichment import cluster_perm as enrich_cluster
from rnamaps.io_rmats import load_rmats_data
from rnamaps.io_vastdb import load_vastdb_data
from rnamaps.logging_utils import log_runtime, setup_logging
from rnamaps.multivalency import plot_multivalency
from rnamaps.permutation import compute_permutation_pvalues
from rnamaps.plots import plot_exon_lengths, plot_heatmap, plot_rna_map
from rnamaps.preprocessing import (
    apply_control_set,
    apply_subsetting,
    autodetect_and_convert_bed_chroms,
    autodetect_and_convert_df_chroms,
    get_ss_bed,
)


def _run_method(method, region_per_exon_df, region_fisher_linegraph,
                exon_categories, region_label, args, rng, smoothing):
    """Dispatch a single enrichment method on one splice-site region.

    Returns
    -------
    EnrichmentResult or None
    """
    if method == 'bootstrap_contrast':
        return enrich_bootstrap.compute(
            region_per_exon_df, exon_categories, region_label,
            rng=rng,
            n_boot=args.n_boot,
            pseudocount=args.pseudocount,
            pseudocount_frac=args.pseudocount_frac,
            bootstrap_control_fixed=args.bootstrap_control_fixed,
        )
    if method == 'cluster_perm':
        return enrich_cluster.compute(
            region_per_exon_df, exon_categories, region_label,
            rng=rng,
            n_perm=args.n_perm,
            cluster_thresh=args.cluster_thresh,
        )
    if method == 'permutation_z':
        plot_df, clusters_df = compute_permutation_pvalues(
            region_per_exon_df, exon_categories, region_label,
            n_perm=args.n_perm, smoothing=smoothing, rng=rng,
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
                 original_counts, window, args, output_dir, FILEname):
    """Render one PDF (or two for bootstrap_contrast) for a method's
    accumulated per-region results."""
    if method == 'bootstrap_contrast':
        subtitle = (f"Bootstrap contrast (B={args.n_boot})"
                    + (' [ctrl fixed]'
                       if args.bootstrap_control_fixed else ''))
        for y_col, ylab in [
            ('delta', 'mean coverage difference (cat - ctrl)'),
            ('log2fc', 'log2 fold change vs control'),
        ]:
            plot_rna_map(
                plot_df, exon_categories, original_counts,
                window, args.all_sites, output_dir, FILEname,
                plot_kind='ribbon', y_col=y_col,
                ci_cols=(f'{y_col}_lo', f'{y_col}_hi'),
                method_name=y_col, ylabel=ylab, subtitle=subtitle,
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
            method_name='cluster_perm',
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
            plot_kind='line', method_name='permutation_z',
        )
        return
    if method == 'fisher':
        plot_rna_map(
            plot_df, exon_categories, original_counts,
            window, args.all_sites, output_dir, FILEname,
            pvalue_method='fisher', plot_kind='line',
            method_name='fisher',
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

        exon_categories = df_rmats.groupby('category').size()
        logging.info("\nExons in each category:")
        logging.info(exon_categories)

        # Validate categories
        if "control" not in exon_categories or exon_categories.loc["control"] == 0:
            logging.error("No control exons found!")
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
                         'FDR', 'dPSI', 'maxPSI',
                         'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']
            save_cols = [c for c in save_cols if c in df_rmats.columns]
            suffix = '_RMATS_with_categories.tsv'
        else:
            save_cols = ['chr', 'exonStart_0base', 'exonEnd', 'strand', 'category',
                         'EVENT', 'GENE',
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

            if getattr(args, 'hg38_chr_autodetect', False):
                xl_bed = autodetect_and_convert_bed_chroms(
                    xl_bed, chroms, args.chr_mapping_file, output_dir)

            middle_3ss = get_coverage_plot(
                xl_bed, middle_3ss_bed, fai, window, exon_categories,
                'middle_3ss', smoothing)
            middle_5ss = get_coverage_plot(
                xl_bed, middle_5ss_bed, fai, window, exon_categories,
                'middle_5ss', smoothing)
            downstream_3ss = get_coverage_plot(
                xl_bed, downstream_3ss_bed, fai, window, exon_categories,
                'downstream_3ss', smoothing)
            upstream_5ss = get_coverage_plot(
                xl_bed, upstream_5ss_bed, fai, window, exon_categories,
                'upstream_5ss', smoothing)

            linegraph_middle_3ss = middle_3ss[0]
            linegraph_middle_5ss = middle_5ss[0]
            linegraph_downstream_3ss = downstream_3ss[0]
            linegraph_upstream_5ss = upstream_5ss[0]

            heatmap_middle_3ss = middle_3ss[1]
            heatmap_middle_5ss = middle_5ss[1]
            heatmap_downstream_3ss = downstream_3ss[1]
            heatmap_upstream_5ss = upstream_5ss[1]

            per_exon_middle_3ss = middle_3ss[2]
            per_exon_middle_5ss = middle_5ss[2]
            per_exon_downstream_3ss = downstream_3ss[2]
            per_exon_upstream_5ss = upstream_5ss[2]

            if not args.all_sites:
                plotting_df = pd.concat([
                    linegraph_upstream_5ss, linegraph_middle_3ss,
                    linegraph_middle_5ss, linegraph_downstream_3ss
                ])
                heat_df = pd.concat([
                    heatmap_upstream_5ss, heatmap_middle_3ss,
                    heatmap_middle_5ss, heatmap_downstream_3ss
                ])
                per_exon_regions = [
                    per_exon_upstream_5ss, per_exon_middle_3ss,
                    per_exon_middle_5ss, per_exon_downstream_3ss,
                ]
            else:
                downstream_5ss = get_coverage_plot(
                    xl_bed, downstream_5ss_bed, fai, window,
                    exon_categories, 'downstream_5ss', smoothing)
                upstream_3ss = get_coverage_plot(
                    xl_bed, upstream_3ss_bed, fai, window,
                    exon_categories, 'upstream_3ss', smoothing)

                linegraph_downstream_5ss = downstream_5ss[0]
                linegraph_upstream_3ss = upstream_3ss[0]
                heatmap_downstream_5ss = downstream_5ss[1]
                heatmap_upstream_3ss = upstream_3ss[1]
                per_exon_downstream_5ss = downstream_5ss[2]
                per_exon_upstream_3ss = upstream_3ss[2]

                plotting_df = pd.concat([
                    linegraph_middle_3ss, linegraph_middle_5ss,
                    linegraph_downstream_3ss, linegraph_downstream_5ss,
                    linegraph_upstream_3ss, linegraph_upstream_5ss
                ])
                heat_df = pd.concat([
                    heatmap_middle_3ss, heatmap_middle_5ss,
                    heatmap_downstream_3ss, heatmap_downstream_5ss,
                    heatmap_upstream_3ss, heatmap_upstream_5ss
                ])
                per_exon_regions = [
                    per_exon_middle_3ss, per_exon_middle_5ss,
                    per_exon_downstream_3ss, per_exon_downstream_5ss,
                    per_exon_upstream_3ss, per_exon_upstream_5ss,
                ]

            # Heatmap (unaffected by enrichment method choice: uses
            # full per-exon data from the original coverage call).
            plot_heatmap(heat_df, exon_categories, window, args.all_sites,
                         output_dir, FILEname)

            # Pair each region's per-exon coverage with its fisher
            # linegraph (used as-is for --enrichment fisher).
            if not args.all_sites:
                fisher_linegraphs = [
                    linegraph_upstream_5ss, linegraph_middle_3ss,
                    linegraph_middle_5ss, linegraph_downstream_3ss,
                ]
            else:
                fisher_linegraphs = [
                    linegraph_middle_3ss, linegraph_middle_5ss,
                    linegraph_downstream_3ss, linegraph_downstream_5ss,
                    linegraph_upstream_3ss, linegraph_upstream_5ss,
                ]

            logging.info("\n" + "=" * 60)
            logging.info(
                f"RUNNING ENRICHMENT METHOD(S): "
                f"{', '.join(enrichment_methods)}"
            )
            logging.info("=" * 60)

            # method -> list of EnrichmentResult per region (in order)
            results_by_method = {m: [] for m in enrichment_methods}

            for region_df, fisher_lg in zip(
                per_exon_regions, fisher_linegraphs
            ):
                if region_df.empty:
                    continue
                region_label = region_df['label'].iloc[0]
                for method in enrichment_methods:
                    res = _run_method(
                        method, region_df, fisher_lg,
                        exon_categories, region_label, args, rng,
                        smoothing,
                    )
                    if res is not None:
                        results_by_method[method].append(res)

            logging.info("\n" + "=" * 60)
            logging.info("PLOTTING RNA MAPS")
            logging.info("=" * 60)

            for method in enrichment_methods:
                method_results = results_by_method[method]
                if not method_results:
                    logging.warning(f"[{method}] No results to plot.")
                    continue

                method_plot_df = pd.concat(
                    [r.plot_df for r in method_results], ignore_index=True
                )
                cluster_frames = [r.clusters_df for r in method_results
                                  if r.clusters_df is not None
                                  and not r.clusters_df.empty]
                method_clusters_df = (pd.concat(cluster_frames,
                                                ignore_index=True)
                                      if cluster_frames else pd.DataFrame())

                method_plot_df.to_csv(
                    f'{output_dir}/{FILEname}_RNAmap_{method}.tsv',
                    sep="\t", index=False,
                )
                if not method_clusters_df.empty:
                    method_clusters_df.to_csv(
                        f'{output_dir}/{FILEname}_RNAmap_{method}_clusters.tsv',
                        sep="\t", index=False,
                    )

                _plot_method(method, method_plot_df, method_clusters_df,
                             exon_categories, original_counts,
                             window, args, output_dir, FILEname)

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

    finally:
        log_runtime(start_time, logger)
        for handler in logger.handlers:
            handler.flush()
        logging.shutdown()

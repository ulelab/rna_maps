"""Top-level RNA map pipeline that orchestrates loading, processing, and plotting."""

import logging
import os
import sys

import numpy as np
import pandas as pd

from rnamaps.coverage import get_coverage_plot
from rnamaps.io_rmats import load_rmats_data
from rnamaps.io_vastdb import load_vastdb_data
from rnamaps.logging_utils import log_runtime, setup_logging
from rnamaps.multivalency import plot_multivalency
from rnamaps.plots import plot_exon_lengths, plot_heatmap, plot_rna_map
from rnamaps.preprocessing import apply_subsetting, get_ss_bed


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
        # Set random seed for reproducibility
        np.random.seed(args.seed)
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

        # Apply subsetting
        if not args.no_subset:
            df_rmats, original_counts = apply_subsetting(
                df_rmats, args.no_constitutive
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

            if not args.all_sites:
                plotting_df = pd.concat([
                    linegraph_upstream_5ss, linegraph_middle_3ss,
                    linegraph_middle_5ss, linegraph_downstream_3ss
                ])
                heat_df = pd.concat([
                    heatmap_upstream_5ss, heatmap_middle_3ss,
                    heatmap_middle_5ss, heatmap_downstream_3ss
                ])
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

            # Save coverage data
            plotting_df.to_csv(
                f'{output_dir}/{FILEname}_RNAmap.tsv', sep="\t", index=False
            )

            # Heatmap
            plot_heatmap(heat_df, exon_categories, window, args.all_sites,
                         output_dir, FILEname)

            # Main RNA map plot
            logging.info("\n" + "=" * 60)
            logging.info("PLOTTING RNA MAPS")
            logging.info("=" * 60)

            plot_rna_map(plotting_df, exon_categories, original_counts,
                         window, args.all_sites, output_dir, FILEname)

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

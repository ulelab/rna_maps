"""Track 2 input loader: VastDB ID lists + EVENT_INFO annotation."""

import logging
import re

import numpy as np
import pandas as pd


def parse_event_info_coordinates(event_info_df):
    """
    Parse genomic coordinates from EVENT_INFO file.
    Returns DataFrame with exon and flanking exon coordinates.
    """
    coords_list = []

    for idx, row in event_info_df.iterrows():
        coord_str = str(row['COORD_o'])
        match = re.match(r'(chr[\w]+):(\d+)-(\d+)', coord_str)
        if not match:
            continue

        ref_co = str(row.get('REF_CO', ''))
        strand_match = re.search(r':([+-])$', ref_co)
        strand = strand_match.group(1) if strand_match else '+'

        upstream_start, upstream_end = np.nan, np.nan
        downstream_start, downstream_end = np.nan, np.nan

        if pd.notna(row.get('CO_C1', np.nan)):
            co_c1_match = re.match(r'chr[\w]+:(\d+)-(\d+)', str(row['CO_C1']))
            if co_c1_match:
                upstream_start = int(co_c1_match.group(1))
                upstream_end = int(co_c1_match.group(2))

        if pd.notna(row.get('CO_C2', np.nan)):
            co_c2_match = re.match(r'chr[\w]+:(\d+)-(\d+)', str(row['CO_C2']))
            if co_c2_match:
                downstream_start = int(co_c2_match.group(1))
                downstream_end = int(co_c2_match.group(2))

        coords_list.append({
            'EVENT': row['EVENT'],
            'GENE': row['GENE'],
            'chr': match.group(1),
            'exonStart_1based': int(match.group(2)),
            'exonEnd': int(match.group(3)),
            'strand': strand,
            'upstreamES': upstream_start,
            'upstreamEE': upstream_end,
            'downstreamES': downstream_start,
            'downstreamEE': downstream_end,
        })

    coords_df = pd.DataFrame(coords_list)

    if len(coords_df) > 0:
        logging.info(f"Parsed coordinates for {len(coords_df)} events")

    return coords_df


def load_vastdb_data(enhanced_file, silenced_file, control_file,
                     constitutive_file, event_info_file, chroms):
    """
    Load VastDB ID lists and look up coordinates from EVENT_INFO.

    Categories are assigned by list membership. Coordinates are converted
    from 1-based (VastDB) to 0-based (BED). CO_C1/CO_C2 are already in
    transcript order, so no strand correction is needed.

    Returns DataFrame with canonical columns for the shared pipeline.
    """
    logging.info("=" * 60)
    logging.info("INPUT MODE: VastDB ID Lists")
    logging.info("=" * 60)

    # Read ID lists
    def read_id_list(filepath, category):
        if filepath is None:
            return []
        with open(filepath) as f:
            ids = [line.strip() for line in f
                   if line.strip() and not line.startswith('#')]
        logging.info(f"Loaded {len(ids)} {category} IDs")
        return ids

    enhanced_ids = read_id_list(enhanced_file, 'enhanced')
    silenced_ids = read_id_list(silenced_file, 'silenced')
    control_ids = read_id_list(control_file, 'control')
    constitutive_ids = read_id_list(constitutive_file, 'constitutive')

    all_ids = enhanced_ids + silenced_ids + control_ids + constitutive_ids
    logging.info(f"Total EVENT IDs (with duplicates): {len(all_ids)}")

    # Check for duplicates
    unique_ids = set(all_ids)
    if len(unique_ids) < len(all_ids):
        n_duplicates = len(all_ids) - len(unique_ids)
        logging.warning(f"\n⚠️  Found {n_duplicates} duplicate EVENT IDs across categories!")

    # Create ID to category mapping (priority: enhanced > silenced > control > constitutive)
    id_to_category = {}
    for eid in constitutive_ids:
        id_to_category[eid] = 'constituitive'

    for eid in control_ids:
        id_to_category[eid] = 'control'

    for eid in silenced_ids:
        id_to_category[eid] = 'silenced'

    for eid in enhanced_ids:
        id_to_category[eid] = 'enhanced'

    logging.info(f"Unique EVENT IDs after deduplication: {len(id_to_category)}")

    # Load EVENT_INFO
    logging.info(f"Loading EVENT_INFO: {event_info_file}")
    event_info = pd.read_csv(event_info_file, sep='\t')
    logging.info(f"EVENT_INFO contains {len(event_info)} events")

    # Parse coordinates
    event_coords = parse_event_info_coordinates(event_info)

    # Filter to our IDs
    event_coords_subset = event_coords[event_coords['EVENT'].isin(all_ids)].copy()

    n_matched = len(event_coords_subset)
    n_total = len(all_ids)
    logging.info(f"Matched {n_matched}/{n_total} IDs to coordinates "
                 f"({100 * n_matched / n_total:.1f}%)")

    # Diagnostic: which IDs were NOT matched
    matched_ids = set(event_coords_subset['EVENT'])
    missing_ids = set(all_ids) - matched_ids

    if len(missing_ids) > 0:
        logging.warning(f"\n⚠️  {len(missing_ids)} IDs not found in EVENT_INFO!")

        missing_by_cat = {}
        for eid in missing_ids:
            cat = id_to_category[eid]
            missing_by_cat[cat] = missing_by_cat.get(cat, 0) + 1

        logging.warning("Missing IDs by category:")
        for cat, count in missing_by_cat.items():
            total_in_cat = len([e for e in all_ids if id_to_category[e] == cat])
            pct = 100 * count / total_in_cat
            logging.warning(f"  {cat}: {count}/{total_in_cat} ({pct:.1f}%)")

        logging.warning("\nFirst 10 missing IDs:")
        for eid in list(missing_ids)[:10]:
            logging.warning(f"  {eid} ({id_to_category[eid]})")

    if n_matched == 0:
        raise ValueError("No EVENT IDs matched to coordinates!")

    # Convert 1-based to 0-based
    event_coords_subset['exonStart_0base'] = event_coords_subset['exonStart_1based'] - 1
    event_coords_subset['upstreamES'] = event_coords_subset['upstreamES'] - 1
    event_coords_subset['downstreamES'] = event_coords_subset['downstreamES'] - 1

    # Assign categories
    event_coords_subset['category'] = event_coords_subset['EVENT'].map(id_to_category)

    # Add FDR placeholder (not used but needed for column contract)
    event_coords_subset['FDR'] = 0.001

    # Filter to valid chromosomes
    event_coords_subset = event_coords_subset[event_coords_subset['chr'].isin(chroms)]

    logging.info(f"Category distribution:\n{event_coords_subset['category'].value_counts()}")

    return event_coords_subset

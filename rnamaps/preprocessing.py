"""Shared preprocessing: subsetting, splice-site BED creation, smoothing."""

import logging
import os

import numpy as np
import pandas as pd


def _load_chrom_mapping(mapping_file):
    if not os.path.exists(mapping_file):
        raise FileNotFoundError(
            f"--hg38_chr_autodetect: chrom mapping file not found: "
            f"{mapping_file}"
        )
    df_map = pd.read_csv(mapping_file, sep='\t', header=None,
                         dtype=str, keep_default_na=False)
    df_map = df_map[(df_map[0] != '') & (df_map[1] != '')]
    ens2gen = dict(zip(df_map[0], df_map[1]))
    gen2ens = dict(zip(df_map[1], df_map[0]))
    return ens2gen, gen2ens


def _pick_chrom_map(src_chroms, target_chroms, ens2gen, gen2ens, label,
                    mapping_file):
    """Pick which mapping direction best maps src_chroms onto target_chroms."""
    overlap = set(src_chroms) & set(target_chroms)
    if len(overlap) > 0 and len(overlap) >= len(set(src_chroms)) * 0.5:
        logging.info(
            f"{label} chrom names already match target "
            f"({len(overlap)}/{len(set(src_chroms))} overlap); no "
            f"conversion needed."
        )
        return None, None

    mapped_ens2gen = {c for c in src_chroms
                      if ens2gen.get(c) in target_chroms}
    mapped_gen2ens = {c for c in src_chroms
                      if gen2ens.get(c) in target_chroms}

    if not mapped_ens2gen and not mapped_gen2ens:
        raise ValueError(
            f"--hg38_chr_autodetect: {label} chrom names do not match "
            f"target and could not be mapped via {mapping_file}. "
            f"{label} chroms (sample): {sorted(set(src_chroms))[:5]}; "
            f"target chroms (sample): {sorted(set(target_chroms))[:5]}"
        )

    if len(mapped_ens2gen) >= len(mapped_gen2ens):
        return ens2gen, "Ensembl -> GENCODE"
    return gen2ens, "GENCODE -> Ensembl"


def autodetect_and_convert_bed_chroms(xl_bed, chroms, mapping_file, output_dir):
    """
    Auto-detect whether the chromosome naming in ``xl_bed`` matches the
    fai chromosomes (``chroms``). If not, convert chrom names using the
    two-column ``mapping_file`` (Ensembl <-> GENCODE) and write a new
    BED file under ``output_dir``. Returns the path to use downstream
    (either the original or the converted file).
    """
    df_bed = pd.read_csv(xl_bed, sep='\t', header=None, comment='#',
                         dtype={0: str})
    bed_chroms = set(df_bed[0].unique())

    ens2gen, gen2ens = _load_chrom_mapping(mapping_file)
    chrom_map, direction = _pick_chrom_map(
        bed_chroms, chroms, ens2gen, gen2ens, "BED", mapping_file)
    if chrom_map is None:
        return xl_bed

    before = len(df_bed)
    df_bed[0] = df_bed[0].map(chrom_map)
    df_bed = df_bed.dropna(subset=[0])
    after = len(df_bed)

    os.makedirs(output_dir, exist_ok=True)
    base = os.path.basename(xl_bed)
    stem, ext = os.path.splitext(base)
    out_path = os.path.join(output_dir, f"{stem}.chrconverted{ext or '.bed'}")
    df_bed.to_csv(out_path, sep='\t', header=False, index=False)

    logging.info(
        f"--hg38_chr_autodetect: converted BED chrom names ({direction}) "
        f"using {mapping_file}; kept {after}/{before} records; "
        f"wrote {out_path}"
    )
    return out_path


def autodetect_and_convert_df_chroms(df, chroms, mapping_file,
                                     chr_col='chr', label='exon'):
    """
    Auto-detect whether ``df[chr_col]`` matches ``chroms`` (fai chroms).
    If not, convert using the Ensembl<->GENCODE ``mapping_file``. Rows
    whose chrom can't be mapped are dropped. Returns the (possibly
    modified) DataFrame.
    """
    src_chroms = set(df[chr_col].astype(str).unique())
    ens2gen, gen2ens = _load_chrom_mapping(mapping_file)
    chrom_map, direction = _pick_chrom_map(
        src_chroms, chroms, ens2gen, gen2ens, label, mapping_file)
    if chrom_map is None:
        return df

    before = len(df)
    df = df.copy()
    df[chr_col] = df[chr_col].astype(str).map(chrom_map)
    df = df.dropna(subset=[chr_col])
    after = len(df)
    logging.info(
        f"--hg38_chr_autodetect: converted {label} chrom names "
        f"({direction}) using {mapping_file}; kept {after}/{before} rows."
    )
    return df


def apply_subsetting(df_rmats, no_constitutive):
    """
    Subset control and constitutive exons to match the largest regulated
    category count. Returns (subsetted df, original_counts dict).
    """
    category_counts = df_rmats['category'].value_counts()
    original_counts = {cat: count for cat, count in category_counts.items()}

    target_count = 0
    if 'enhanced' in category_counts and 'silenced' in category_counts:
        target_count = max(category_counts['enhanced'], category_counts['silenced'])
    elif 'enhanced' in category_counts:
        target_count = category_counts['enhanced']
    elif 'silenced' in category_counts:
        target_count = category_counts['silenced']

    # Subset control
    if 'control' in category_counts and category_counts['control'] > target_count > 0:
        control_indices = df_rmats[df_rmats['category'] == 'control'].index
        control_keep = np.random.choice(control_indices, target_count, replace=False)
        drop_mask = df_rmats.index.isin(control_indices) & ~df_rmats.index.isin(control_keep)
        df_rmats = df_rmats[~drop_mask]
        logging.info(f"Subsetted control exons from {category_counts['control']} to {target_count}")

    # Subset constitutive
    if (not no_constitutive
            and 'constitutive' in category_counts
            and category_counts['constitutive'] > target_count > 0):
        const_indices = df_rmats[df_rmats['category'] == 'constitutive'].index
        const_keep = np.random.choice(const_indices, target_count, replace=False)
        drop_mask = df_rmats.index.isin(const_indices) & ~df_rmats.index.isin(const_keep)
        df_rmats = df_rmats[~drop_mask]
        logging.info(f"Subsetted constitutive exons from "
                     f"{category_counts['constitutive']} to {target_count}")

    return df_rmats, original_counts


def get_ss_bed(df, pos_col, neg_col):
    """
    Create BED file for splice sites (handles strand orientation).

    pos_col: column to use for + strand (one end of the exon)
    neg_col: column to use for - strand (other end of same exon)

    Both columns should reference the SAME exon. The upstream/downstream
    strand correction happens at load time, not here.
    """
    df = df.copy()

    df['exon_id'] = (
        df['category'] + "_" +
        df['chr'].astype(str) + ":" +
        df['exonStart_0base'].astype(str) + "-" +
        df['exonEnd'].astype(str) + ";" +
        df['strand'].astype(str)
    )

    ss_pos = df.loc[df['strand'] == "+",
                    ['chr', pos_col, pos_col, 'exon_id', 'FDR', 'strand']]
    ss_pos.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_pos.start = ss_pos.start.transform(lambda x: x - 1)

    ss_n = df.loc[df['strand'] == "-",
                  ['chr', neg_col, neg_col, 'exon_id', 'FDR', 'strand']]
    ss_n.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
    ss_n.end = ss_n.end.transform(lambda x: x + 1)

    ss = pd.concat([ss_pos, ss_n])

    return ss


def smooth_coverage(df, window_size=10, std=2):
    """Smooth coverage data using rolling Gaussian window."""
    result = df.copy()
    groups = []

    for (exon_id, label), group_df in result.groupby(['exon_id', 'label']):
        if len(group_df) < 3:
            groups.append(group_df)
            continue

        group_sorted = group_df.sort_values('position')

        if len(group_sorted) >= window_size:
            values = group_sorted['coverage'].values
            s = pd.Series(values)

            smoothed = s.rolling(
                window=window_size,
                center=True,
                win_type='gaussian'
            ).mean(std=std)

            smoothed = smoothed.fillna(s)
            group_sorted['coverage'] = smoothed.values

        groups.append(group_sorted)

    return pd.concat(groups)

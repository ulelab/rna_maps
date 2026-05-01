## RNA maps
Authors: charlotte.capitanchik@crick.ac.uk; leomwilkinson@gmail.com; aram.amalietti@gmail.com

---

## Overview

`rnamaps` generates RNA maps showing positional enrichment of RBP binding (from CLIP data) around regulated exons. It supports two input modes:

1. **rMATS mode** — takes rMATS differential splicing output and auto-categorises exons from dPSI/FDR thresholds
2. **VastDB mode** — takes pre-curated VastDB EVENT ID lists with categories already assigned

Both modes feed into the same analysis pipeline: splice site BED creation, CLIP coverage calculation, label-permutation test for positional enrichment, RNA map plotting, per-exon heatmaps, and exon length distributions.

---

## Quick Start

Create a conda environment and install the package:

```
conda env create -f environment.yml
conda activate rnamaps
pip install -e .
```

After installation, the `rnamaps` command is available. You can also run the package directly with `python -m rnamaps`.

Small test for rMATS mode:
```
python rna_maps.py \
-i test/chr21_PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt \
-x test/chr21_hela_ptbp1_iclip_sorted_merged.bed \
-f test/homosapien-hg37-chr21.fa \
-fi test/homosapien-hg37-chr21.fa.fai
```
Small test for VastDB mode:
```
rnamaps \
  --vastdb_mode \
  --vastdb_enhanced test/vast-tools/enhanced_ids.chr21.txt \
  --vastdb_silenced test/vast-tools/silenced_ids.chr21.txt \
  --vastdb_control test/vast-tools/control_ids.chr21.txt \
  --vastdb_constitutive test/vast-tools/constitutive_ids.chr21.txt \
  --vastdb_annotation test/vast-tools/EVENT_INFO-hg38.chr21.tab \
  -x test/chr21_hela_ptbp1_iclip_sorted_merged.bed \
  -f test/homosapien-hg37-chr21.fa \
  -fi test/homosapien-hg37-chr21.fa.fai \
  -o output/ \
  -p PTBP1
```

### rMATS mode

```
rnamaps \
  -i SE.MATS.JCEC.txt \
  -x CLIP_crosslinks.bed \
  -f genome.fa \
  -fi genome.fa.fai \
  -o output/ \
  -p PTBP1
```

### VastDB mode

```
rnamaps \
  --vastdb_mode \
  --vastdb_enhanced enhanced_ids.txt \
  --vastdb_silenced silenced_ids.txt \
  --vastdb_control control_ids.txt \
  --vastdb_constitutive constitutive_ids.txt \
  --vastdb_annotation EVENT_INFO-hg38.tab \
  -x CLIP_crosslinks.bed \
  -f genome.fa \
  -fi genome.fa.fai \
  -o output/ \
  -p PTBP1
```

---

## Input modes

### rMATS mode (`-i`)

Accepts rMATS quantified files for cassette exons (e.g. `SE.MATS.JCEC.txt`). Categories are assigned automatically based on dPSI and FDR thresholds (see [Definitions](#definitions) below). The script deduplicates exons involved in multiple events by keeping the most extreme dPSI value per exon.

If your condition is RBP knockdown, run your comparison as condition − control, such that definitions of enhanced and silenced are correct. If your condition is RBP overexpression, run as control − condition. In the generic comparison group1 − group2, "enhanced" and "silenced" are defined in reference to group2 relative to group1.

**Minus-strand correction:** rMATS labels flanking exons as "upstream" and "downstream" by genomic coordinate (lower = upstream), which is inverted for minus-strand genes. The script automatically swaps these at load time so that upstream/downstream always refer to transcript order throughout the pipeline.

### VastDB mode (`--vastdb_mode`)

Accepts four plain text files of VastDB EVENT IDs (one per line, `#` comments allowed), each representing a pre-assigned category. Genomic coordinates are looked up from the VastDB `EVENT_INFO` annotation file using the `COORD_o`, `CO_C1`, and `CO_C2` columns.

This mode is useful when categories come from VAST-TOOLS `vast diff` output (or any other source of splicing quantification), where you have already applied your own thresholds to define enhanced, silenced, control, and constitutive exons. No rMATS dependency is required.

| Feature | rMATS mode | VastDB mode |
|---|---|---|
| Splicing input | rMATS `SE.MATS.JCEC.txt` | VastDB EVENT ID lists |
| Category assignment | Computed from dPSI / FDR thresholds | Pre-assigned by user |
| Coordinate source | rMATS columns | VastDB `EVENT_INFO` file |
| Coordinate system | 0-based (BED) | 1-based (converted automatically) |
| rMATS dependency | Required | Not required |
| Expression filtering | From junction counts | Not applicable |

---

## Preparing inputs

### CLIP crosslinks

iCLIP or eCLIP crosslink sites in BED format. Single-nucleotide resolution crosslink files (`.genome.xl.bed`) are preferred over peak files for RNA maps.

### Genome reference

Genome FASTA (`.fa`) and FASTA index (`.fa.fai`) files. The genome build must match both the splicing data and the CLIP data — e.g. hg19 for older rMATS datasets, hg38 for VastDB.

### rMATS input

Standard rMATS output for skipped exons, e.g. `SE.MATS.JCEC.txt`. Must contain columns: `chr`, `exonStart_0base`, `exonEnd`, `strand`, `FDR`, `IncLevelDifference`, `IncLevel1`, `IncLevel2`, `upstreamES`, `upstreamEE`, `downstreamES`, `downstreamEE`.

### VastDB inputs

- **EVENT ID lists**: Plain text files with one VastDB EVENT ID per line (e.g. `HsaEX0012345`). Lines beginning with `#` are ignored. At least one list file must be provided; typically all four categories are supplied.
- **VastDB annotation file**: The `EVENT_INFO-*.tab` file from VastDB. Must contain `EVENT`, `GENE`, `COORD_o`, `REF_CO`, `CO_C1`, and `CO_C2` columns. The version in your `event_lists/` directory (with the full `CO_C1`/`CO_C2` columns) is required — not the minimal version.

---

## Usage

```
rnamaps -h

Input modes (mutually exclusive, one required):
  -i, --inputsplice     rMATS differential splicing file (rMATS mode)
  --vastdb_mode         Use VastDB ID lists mode

VastDB mode options:
  --vastdb_enhanced     Enhanced exon IDs (one per line)
  --vastdb_silenced     Silenced exon IDs (one per line)
  --vastdb_control      Control exon IDs (one per line)
  --vastdb_constitutive Constitutive exon IDs (one per line)
  --vastdb_annotation   VastDB EVENT_INFO file

Required arguments (both modes):
  -x, --inputxlsites    CLIP crosslinks in BED file format
  -f, --genomefasta     Genome FASTA file (.fa)
  -fi, --fastaindex     Genome FASTA index (.fai)

Optional arguments:
  -o, --outputpath      Output folder [DEFAULT: current directory]
  -w, --window          Window around splice sites [DEFAULT: 300]
  -s, --smoothing       Smoothing window [DEFAULT: 15]
  --seed                Random seed for reproducible permutations / subsetting [DEFAULT: 42]
  -nc, --no_constitutive  Exclude constitutive category
  -ns, --no_subset      Disable subsetting of control/constitutive exons
                        (subsetting is auto-disabled when --permute is on)
  -ao, --all_sites      Include all 6 splice sites (default: 4 core sites)
  -p, --prefix          Prefix for output files

Permutation test options:
  --permute / --no-permute
                        Use label-permutation test for p-values [DEFAULT: --permute].
                        --no-permute switches to a per-position Fisher's exact test
                        and re-enables exon subsetting.
  --n_perm              Number of label permutations [DEFAULT: 1000]

rMATS mode thresholds:
  -mc, --minctrl        Minimum dPSI for control events [DEFAULT: -0.05]
  -xc, --maxctrl        Maximum dPSI for control events [DEFAULT: 0.05]
  -xi, --maxincl        Maximum PSI for control (above = constitutive) [DEFAULT: 0.9]
  -xf, --maxfdr         Maximum FDR for regulated events [DEFAULT: 0.1]
  -xe, --maxenh         Maximum dPSI for enhanced exons [DEFAULT: -0.05]
  -ms, --minsil         Minimum dPSI for silenced exons [DEFAULT: 0.05]

Multivalency analysis:
  -v, --multivalency    Run multivalency analysis (requires germs.R)
  -g, --germsdir        Directory containing germs.R [DEFAULT: current directory]
```

---

## Outputs

Both modes produce the same set of output files:

| File | Description |
|---|---|
| `{prefix}_RMATS_with_categories.tsv` or `{prefix}_VastDB_with_categories.tsv` | Categorised exons with coordinates |
| `{prefix}_RNAmap.tsv` | Per-position coverage and enrichment data |
| `{prefix}_RNAmap_-log10pvalue.pdf` | RNA map line plot |
| `{prefix}_heatmap.pdf` | Per-exon binary coverage heatmap |
| `{prefix}_totalExonsCovered.tsv` | Count of exons with CLIP signal per region and category |
| `{prefix}_exon_length.pdf` | Exon length distributions by category |
| `execution_*.log` | Run log with timing and category counts |

With `--multivalency` (rMATS mode):

| File | Description |
|---|---|
| `{prefix}_RNAmap_multivalency.tsv` | Multivalency scores per position |
| `{prefix}_RNAmap_multivalency.pdf` | Multivalency plot |
| `{prefix}_RNAmap_TOP10KMER_multivalency.tsv` | Top kmer multivalency scores |
| `{prefix}_RNAmap_silencedKMER_multivalency.pdf` | Silenced exon kmer plot |
| `{prefix}_RNAmap_enhancedKMER_multivalency.pdf` | Enhanced exon kmer plot |

---

## Multivalency analysis

Multivalency analysis adds run time and requires the Ule lab's GeRMs package. It is optional and enabled with the `-v` flag.

To install GeRMs, clone the repository and run from within it (requires R devtools):

```
R -e 'devtools::install()'
```

GeRMs requires: Biostrings, parallel, logger, and optparse.

When running RNA maps with multivalency, provide the location of the germs repo with `-g`:

```
rnamaps \
  -i SE.MATS.JCEC.txt \
  -x CLIP_crosslinks.bed \
  -f genome.fa \
  -fi genome.fa.fai \
  -v -g ../germs
```

To create a multivalency map without CLIP data, run the above command without the `-x` flag.

---

## Definitions

### Event types

**Control**: An event that doesn't change in inclusion (PSI) in this RBP knockdown, but might in another circumstance.

```
dPSI   ( -1 <---------- - 0.05xxxxx0xxxxx0.05----------> 1 )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9------> 1 )
FDR    (  0 xxxxx0.1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 1 )
```

**Constitutive**: An event that doesn't change in inclusion in this knockdown and is unlikely to change in another circumstance. Defined as a control event with maximum inclusion (PSI) > 0.9.

```
dPSI   ( -1 <---------- - 0.05xxxxx0xxxxx0.05----------> 1 )
maxPSI (  0 ----------------------------------0.9xxxxxx> 1 )
FDR    (  0 xxxxx0.1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 1 )
```

**Enhanced**: An event that is *less* included in RBP knockdown, suggesting the RBP *promotes/enhances* inclusion of the event.

```
dPSI   ( -1 <xxxxxxxxxxx- 0.05-----0-----0.05----------> 1 )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9xxxxxx> 1 )
FDR    (  0 xxxxx0.1-----------------------------------> 1 )
```

**Silenced**: An event that is *more* included in RBP knockdown, suggesting the RBP *represses/silences* inclusion of the event.

```
dPSI   ( -1 <---------- - 0.05-----0-----0.05xxxxxxxxxx> 1 )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9xxxxxx> 1 )
FDR    (  0 xxxxx0.1-----------------------------------> 1 )
```

**Enhanced/Silenced rest** (rMATS mode only): A silenced or enhanced event where the FDR does not fall below the threshold.

```
dPSI   (            As in silenced or enhanced             )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9xxxxxx> 1 )
FDR    (  0 -----0.1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 1 )
```

### Hierarchy (rMATS mode)

An exon may be involved in multiple events. To avoid plotting it multiple times, a hierarchy is applied:

1. If an exon meets criteria for silenced or enhanced, this is designated. If criteria for both are met, the most extreme dPSI value is preferred.
2. Of remaining exons, if they meet criteria for enhanced/silenced rest, this is designated.
3. Of remaining exons, if they meet criteria for constitutive, this is designated.
4. Of remaining exons, if they meet criteria for control, this is designated.

### VastDB mode categories

In VastDB mode, categories are assigned by the user before running the script. Typical thresholds when using VAST-TOOLS `vast diff` output:

| Category | E[dPsi] | MV[dPsi]_at_0.95 | maxPSI |
|---|---|---|---|
| Enhanced | < −0.10 | > 0.05 | — |
| Silenced | > 0.10 | > 0.05 | — |
| Control | \|E[dPsi]\| < 0.05 | — | < 0.9 |
| Constitutive | \|E[dPsi]\| < 0.05 | — | ≥ 0.9 |

---

## Dependencies

These are the versions the script was developed with (pandas >= 1 introduced breaking changes):

```
python=3.11
pandas
numpy
pybedtools
matplotlib
seaborn
scipy
```

All dependencies are specified in `environment.yml` (for conda) and `pyproject.toml` (for pip).

---

## Reproducibility

The `--seed` flag (default: 42) controls the random seed used for both the
label-permutation test and (when `--no-permute` is in effect) the random
subsetting of control / constitutive exons. Setting the same seed produces
identical results across runs.

---

## Statistical test

By default (`--permute`) `rnamaps` calls position-wise enrichment with a
**label-permutation test**. A per-position Fisher's exact test is also
available via `--no-permute`.

### What gets tested

For each splice-site region (e.g. `middle_3ss`) and each non-control
category `c` ∈ {enhanced, silenced, constitutive}, the pipeline:

1. Pools the `n_c` exons in `c` with the `n_ctrl` control exons.
2. Computes the observed test statistic per position `p`:

   `T_obs(p) = mean_coverage_c(p) − mean_coverage_ctrl(p)`

   This is a signed effect size (positive ⇒ enriched relative to control)
   that is well-defined even when the control coverage is zero.
3. Repeats `B = --n_perm` times: randomly relabels `n_c` of the
   `n_c + n_ctrl` pooled exons as "category", recomputes
   `T_perm(p)`, and tallies how often `|T_perm(p)| ≥ |T_obs(p)|`.
4. Returns a **two-sided empirical p-value** per position, with the
   standard `(1 + #ge) / (B + 1)` correction so p is bounded away from 0.
5. Plots **signed −log10(p)**, with the sign taken from the direction of
   `T_obs(p)`. The signed −log10(p) curve is then Gaussian-smoothed
   (controlled by `-s, --smoothing`) for visual clarity.

### Why permutation is the default

- **Imbalanced n is fine.** The permutation null is constructed from the
  actual data, so unequal numbers of regulated vs. control exons do not
  inflate the test statistic. As a result, control / constitutive
  subsetting is automatically disabled
  when `--permute` is on — all controls are used.
- **Distributional assumptions are weaker.** Fisher's exact test treats
  each position independently as a 2×2 contingency table; permutation
  only assumes exchangeability of exon labels under the null.
- **Reproducible:** the random number generator is seeded by `--seed`.

Under `--no-permute` the pipeline instead computes a Fisher's exact test
per position against the control set; in that mode control and
constitutive exons are subset to match the largest regulated category
(seeded by `--seed`).

### Choosing `--n_perm`

`--n_perm` (B) sets the **resolution of the empirical p-value**, not the
statistical power of the test. The smallest p you can resolve is
`1 / (B + 1)`. Reasonable choices:

| `--n_perm` | Smallest reportable p | Typical use |
|---|---|---|
| 100 | ≈ 0.01 | quick sanity / smoke test |
| **1000** (default) | ≈ 0.001 | normal usage |
| 10000 | ≈ 0.0001 | only if you want to plot very small p |

**The number of regulated exons does *not* dictate `B`.** It determines
the *power* of the test (whether real signal crosses your significance
threshold), which raising `B` cannot fix. The only edge case where exon
counts limit `B` is when the total number of unique relabellings,
`C(n_c + n_ctrl, n_c)`, is smaller than `B` — for realistic RNA-map sizes
(dozens of regulated exons, ≥50 controls) this combinatorial space is
enormous, so `B` and exon count are effectively independent.

### Caveats

- Smoothing is applied to the signed −log10(p) curve. This is a cosmetic step; the
  underlying p-values in `{prefix}_RNAmap.tsv` are unsmoothed.
- The test does not control for any covariates (exon length, GC content,
  expression). If you need covariate matching, pre-filter your control
  set accordingly.
- No multiple-testing correction is applied across the hundreds of
  correlated positions in each region. Treat individual peaks as
  exploratory; biological replication of the curve shape is the strongest
  evidence.

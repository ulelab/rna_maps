## RNA maps
Authors: charlotte.capitanchik@crick.ac.uk; leomwilkinson@gmail.com; aram.amalietti@gmail.com

---

## Overview

`rnamaps` generates RNA maps showing positional enrichment of RBP binding (from CLIP data) around regulated exons. It supports two input modes:

1. **rMATS mode** — takes rMATS differential splicing output and auto-categorises exons from dPSI/FDR thresholds
2. **VastDB mode** — takes pre-curated VastDB EVENT ID lists with categories already assigned

Both modes feed into the same analysis pipeline: splice site BED creation, CLIP coverage calculation, optional per-exon binarisation (each exon × position cell becomes 0/1 = "has at least one crosslink at this base" — applied by default for `bootstrap_contrast`, `cluster_perm`, `permutation_z`, `fisher`; **not** applied by `roc_auc`), one or more enrichment analyses (bootstrap contrast / cluster permutation / legacy Fisher / legacy permutation z-score / ROC-AUC), RNA map plotting, per-exon heatmaps, and exon length distributions.

The `roc_auc` method additionally treats the BED score column as a continuous predictor (configurable via `--xl_score`) — useful for AI prediction tracks or any per-base score, not just integer crosslink counts.

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
rnamaps \
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
| Expression matching (`--gene_tpm`) | Supported | Supported |

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
  -s, --smoothing       Smoothing window [DEFAULT: 15]. Centred
                        Gaussian rolling mean applied to each bootstrap
                        iteration's per-position delta/log2fc (and to
                        permutation_z / fisher -log10(p)).
  --seed                Random seed for reproducible permutations / subsetting [DEFAULT: 42]
  -nc, --no_constitutive  Exclude constitutive category
  -ns, --no_subset      Disable subsetting of control/constitutive exons
                        (subsetting is auto-disabled for any non-fisher method)
  -ao, --all_sites      Include all 6 splice sites (default: 4 core sites)
  -p, --prefix          Prefix for output files
  --enrichment          One or more of {bootstrap_contrast, cluster_perm,
                        permutation_z, fisher, roc_auc}
                        [DEFAULT: bootstrap_contrast].
                        See "Enrichment methods" below.
  --y_axis              Y-axis for legacy permutation_z plot:
                        log10p (default) or zscore.

Control-set hygiene options:
  --control_set         {default, strict, constitutive_only} [DEFAULT: default].
                        - strict: tighten control to |dPSI|<--control_max_dpsi
                          AND FDR>--control_min_fdr.
                        - constitutive_only: drop the original control
                          category and relabel constitutive -> control.
  --control_max_dpsi    Strict mode: max |dPSI| [DEFAULT: 0.01]
  --control_min_fdr     Strict mode: min FDR    [DEFAULT: 0.5]

Expression-matching options:
  --gene_tpm            Optional 2-column table (gene_id, tpm). When set,
                        control/constitutive pools are expression-matched
                        to regulated exons.
  --tpm_n_bins          Quantile bins for regulated log-TPM [DEFAULT: 10]
  --tpm_pseudocount     Pseudocount for log10(TPM+p) [DEFAULT: 1.0]
  --tpm_min_tpm         Drop rows with TPM below this value [DEFAULT: 0.0]
  --no_match_constitutive
                        Match only control (leave constitutive unchanged)

Permutation test options:
  --permute / --no-permute
                        Legacy switch retained for backward compatibility.
                        --no-permute is equivalent to --enrichment fisher
                        when --enrichment is not given.
  --n_perm              Number of label permutations [DEFAULT: 1000]
                        (used by permutation_z and cluster_perm)

Bootstrap contrast options (--enrichment bootstrap_contrast):
  --n_boot              Bootstrap iterations [DEFAULT: 1000]
  --bootstrap_control_fixed
                        Treat control mean as a constant (skip resampling
                        control). Equivalent to full bootstrap up to
                        negligible variance when n_ctrl >> n_c.
  --shrinkage           {magnitude, pseudocount, none}
                        How log2fc is regularised so sparse positions
                        don't blow up [DEFAULT: magnitude].
                        See "Shrinkage for log2fc" below for a plain-
                        language explanation.
  --shrinkage_scale     Only used with --shrinkage magnitude. Rate
                        scale tau at which shrinkage transitions from
                        heavy to light [DEFAULT: 0.05, i.e. 5% rate].
  --pseudocount         Only used with --shrinkage pseudocount. Override
                        adaptive log2FC pseudocount with a fixed value
                        [DEFAULT: adaptive].
  --pseudocount_frac    Only used with --shrinkage pseudocount. Adaptive
                        pseudocount fraction of the regional control
                        coverage median [DEFAULT: 0.01].

Cluster-permutation options (--enrichment cluster_perm):
  --cluster_thresh      Cluster-defining |t| threshold [DEFAULT: 2.0]

BED score handling and ROC/AUC options:
  --xl_score            One or more of {ignore, raw, per_transcript_zscore}.
                        How to use BED column 5 of -x.
                        Pass several to sweep modes in one invocation.
                        [DEFAULT: ignore]. See "ROC / AUC analysis" below.
  --roc_aggregator      {mean, max, both}. How per-exon signal is collapsed
                        over the window for the per-region ROC curve.
                        [DEFAULT: both]
  --roc_n_perm          Number of label permutations for the AUC null.
                        0 = skip; >0 emits per-position and per-region
                        permutation p-values. [DEFAULT: 0]

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

Both modes produce a common set of files plus per-method enrichment outputs.

Always produced:

| File | Description |
|---|---|
| `{prefix}_RMATS_with_categories.tsv` or `{prefix}_VastDB_with_categories.tsv` | Categorised exons with coordinates |
| `{prefix}_heatmap.pdf` | Per-exon binary coverage heatmap |
| `{prefix}_totalExonsCovered.tsv` | Count of exons with CLIP signal per region and category |
| `{prefix}_exon_length.pdf` | Exon length distributions by category |
| `execution_*.log` | Run log with timing and category counts |

Per enrichment method (one set per method passed to `--enrichment`):

| Method | Files |
|---|---|
| `bootstrap_contrast` | `{prefix}_RNAmap_bootstrap_contrast.tsv`, `{prefix}_RNAmap_delta.pdf`, `{prefix}_RNAmap_log2fc.pdf` |
| `cluster_perm` | `{prefix}_RNAmap_cluster_perm.tsv`, `{prefix}_RNAmap_cluster_perm_clusters.tsv`, `{prefix}_RNAmap_cluster_perm.pdf` |
| `permutation_z` | `{prefix}_RNAmap_permutation_z.tsv`, `{prefix}_RNAmap_permutation_z.pdf` |
| `fisher` | `{prefix}_RNAmap_fisher.tsv`, `{prefix}_RNAmap_fisher.pdf` |
| `roc_auc` | `{prefix}_RNAmap_roc_auc.tsv` (per-position AUC table), `{prefix}_RNAmap_roc_auc.pdf` (per-position AUC line), `{prefix}_RNAmap_roc_auc_region_auc.tsv` (headline per-region AUCs), `{prefix}_RNAmap_roc_auc_roc_curves.tsv` (FPR/TPR points), `{prefix}_RNAmap_roc_auc_curves_mean.pdf` and `{prefix}_RNAmap_roc_auc_curves_max.pdf` (per-region ROC plots, one per aggregator) |

When multiple `--xl_score` modes are passed in one invocation, each per-method file above gets an `_xlscore-{mode}` suffix (e.g. `{prefix}_RNAmap_roc_auc_xlscore-raw.pdf`). The xl_score-invariant outputs (heatmap, exon length, totalExonsCovered, RMATS_with_categories) are still written once.

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

The `--seed` flag (default: 42) controls the random seed used by all
methods that resample (bootstrap_contrast, cluster_perm, permutation_z)
and by the legacy random subsetting of control / constitutive exons.
Setting the same seed produces identical results across runs.

---

## Enrichment methods

`rnamaps` supports several enrichment methods, selected (one or more) via
`--enrichment`. The default is `bootstrap_contrast`.

### Decision table

| Method | Output | Use when |
|---|---|---|
| `bootstrap_contrast` *(default)* | per-position `delta` and `log2fc` with bootstrap CI ribbons | You want signed effect size + uncertainty, not a p-value. Both contrasts are `n_c`-invariant, so they reward higher per-exon coverage rate rather than larger categories. Especially good when `n_ctrl >> n_c` or you want to compare maps from different libraries. |
| `cluster_perm` | per-cluster p-values, plot of Welch t with significance bars | You want significance per peak with FWER control across positions and don't want a saturating p-value floor. |
| `permutation_z` | signed -log10(p) and z-score per position | Legacy: per-position significance via label-permutation z-score. Useful for comparing to historical results. |
| `fisher` | signed -log10(p) per position | Legacy: per-position Fisher's exact test. Subsetting of control/constitutive is auto-enabled (others auto-disable it). |
| `roc_auc` | per-position AUC line + per-region ROC curves and headline AUCs | You're evaluating a continuous predictor (AI binding probability, peak height, signal density) rather than a binary "is this exon bound" indicator. Rank-based, so robust to score distribution; particularly good with `--xl_score raw` / `per_transcript_*` to actually use BED column 5 instead of treating it as 1. |

You can pass multiple methods at once: `--enrichment bootstrap_contrast cluster_perm` produces the union of their outputs.

### `bootstrap_contrast` (default)

**Per-exon binarisation.** Before any per-exon enrichment method runs,
the per-exon coverage matrix is binarised to 0/1 per (exon, base): each
cell answers *"did this exon have at least one crosslink at this base?"*
rather than reporting a raw read count. This makes
`mean_cov_c(p)` the **fraction of category exons positive at position
`p`**, and prevents a handful of highly-expressed exons from dominating
per-exon statistics through their read counts alone. The legacy
`fisher` path still uses raw read totals; only the per-exon methods
(`bootstrap_contrast`, `cluster_perm`, `permutation_z`) operate on the
binarised matrix.

For each non-control category `c` vs. control, let `n_c` be the number
of exons in `c` and `n_ctrl` the number of control exons:

1. Resample `n_c` exons with replacement from `c` and `n_ctrl` from
   control, `B = --n_boot` times.
2. Per iteration `b`, compute `mean_cov_c^b(p)` and `mean_cov_ctrl^b(p)`
   (a fraction in `[0, 1]` after binarisation: "fraction of exons in
   the resample positive at position `p`").
3. Report two contrasts per position with their `(2.5, 97.5)`
   percentile bands across `b`:

   - `delta(p) = mean_cov_c(p) - mean_cov_ctrl(p)` — additive, rate
     scale. Invariant to `n_c`: it rewards a higher *fraction* of
     category exons being positive at `p`, not larger categories.
   - `log2fc(p) = log2((mean_cov_c(p) + ε) / (mean_cov_ctrl(p) + ε))` —
     multiplicative scale, library-size invariant.

   Where the CI band excludes 0, the contrast is "significant" in a
   bootstrap sense, with no p-value involved.

**Smoothing.** When `--smoothing > 1` (default 15), each bootstrap
iteration's per-position `delta_b` and `log2fc_b` are convolved with a
centred Gaussian-weighted rolling mean **before** taking the across-
iteration mean and percentiles. Smoothing first makes the CI band the
correct uncertainty band for the smoothed estimator; smoothing the
percentiles afterwards would understate uncertainty at sharp features.

**Class imbalance** (e.g. `n_ctrl = 10000`, `n_c = 500`) is handled
cleanly: the bootstrap variance of `mean_cov_ctrl(p)` is negligible, so
the CI on the contrast reflects the category-side variance, which is the
right behaviour. For a 5-10× speed-up at large `n_ctrl`, pass
`--bootstrap_control_fixed` to skip resampling control entirely (treat
its mean as a constant). The pipeline auto-suggests this when
`n_ctrl ≥ 20 × n_c`.

**Shrinkage for log2fc** is controlled by `--shrinkage`. The default is
`beta_binomial`. See "Shrinkage for log2fc" below for a plain-language
explanation of what each mode does and when to use it.

**Important**: the bootstrap reflects sampling uncertainty in the
estimator. It does *not* correct for *bias* from contaminated controls
(silently regulated exons hiding in the control pool). Use
`--control_set` for that.

### `cluster_perm` (Maris-Oostenveld cluster-mass permutation)

For each non-control category vs. control:

1. Per position compute Welch t-statistic
   `t_obs(p) = (mean_c(p) - mean_ctrl(p)) / sqrt(var_c(p)/n_c + var_ctrl(p)/n_ctrl)`.
2. Find contiguous runs where `|t_obs(p)| > --cluster_thresh` (default 2.0).
   Cluster mass is `sum(t_obs)` over the run.
3. Repeat `B = --n_perm` times: shuffle category/control labels,
   recompute `t_perm`, find clusters, record `max(|cluster_mass|)`.
4. Cluster p:
   `p_cluster = (1 + #(null_max ≥ |obs_mass|)) / (B + 1)`.

This controls family-wise error across positions while respecting the
position-to-position correlation that pointwise tests ignore. Output
includes a per-position table (`t_obs`, `in_sig_cluster`,
`cluster_pvalue`) and a cluster summary (`start_pos, end_pos, mass,
cluster_pvalue, sign`). Plot shows the `t_obs` curve plus horizontal
bars beneath the axis at significant clusters (p ≤ 0.05).

### `permutation_z` (legacy z-score against permutation null)

The previous default test. For each position the observed
`T_obs(p) = mean_c(p) - mean_ctrl(p)` is standardised against the
permutation null mean and standard deviation, and a two-sided p-value
is reported via a normal-tail approximation. The reported `-log10(p)`
is unbounded (does not saturate at `1/(B+1)`).

### `fisher` (legacy per-position Fisher's exact test)

For each position, a 2×2 contingency table of {covered, not covered} ×
{category, control} is tested with Fisher's exact. The signed
`-log10(p)` is reported, with sign taken from the per-position
fold-change. With `--enrichment fisher` (or legacy `--no-permute`),
control and constitutive exons are randomly subset to match the largest
regulated category (seeded by `--seed`).

---

## ROC / AUC analysis (`--enrichment roc_auc`)

The other enrichment methods ask "does CLIP signal differ between
regulated and control exons at this position?". `roc_auc` asks the
reciprocal question: **"how good is the CLIP / AI score as a *predictor*
of whether an exon is regulated?"** That's a binary-classification
framing and the natural summary is an ROC curve with its AUC.

Two complementary outputs are produced. They use the same per-exon
coverage matrix as every other method, but interpret it differently.

### Setup: what is the "score" and what is the "label"?

After the coverage step you have a per-exon × per-position matrix
`X` of shape `(n_exons, 2*window + 1)`. The entries of `X` are
controlled by `--xl_score`:

- `--xl_score ignore` (default, legacy): `X[i, p]` is the number of
  crosslinks from the BED that overlap exon `i` at base `p`. Same
  matrix the other methods see; usually 0 or 1 for typical iCLIP
  inputs, larger for heavily-sequenced eCLIP.
- `--xl_score raw`: `X[i, p]` is the **sum of BED column 5 values**
  for all crosslinks overlapping exon `i` at base `p`. For an iCLIP
  BED where column 5 is the merged crosslink count, this gives you
  read depth at the base; for an AI track where column 5 is a
  per-base probability, this gives you the AI score at the base.
- `--xl_score per_transcript_zscore`: same as `raw` but with the BED
  entries' scores **renormalised within each BED column 4 (`name`)
  group** (z-scored per transcript) before they're accumulated onto
  the exon matrix. See "How `--xl_score` impacts these calculations"
  below.

For each non-control category `c` (e.g. `enhanced`) vs. control:

- **Positive class** = exons in category `c` → label `y_i = 1`.
- **Negative class** = exons in control → label `y_i = 0`.
- For the per-position AUC at position `p`, the **score** is the row
  of `X` at that single column: `s_i = X[i, p]`.
- For the per-region ROC, the score is an aggregate over all
  positions: `s_i = mean(X[i, :])` or `s_i = max(X[i, :])` (selected
  with `--roc_aggregator`).

`roc_auc` does **not** binarise. The matrix entries flow into the
rank statistic as-is — that's how `--xl_score raw` differs from
`--xl_score ignore`.

### How AUC is calculated per position

For a single position `p` we have a vector of `n_cat + n_ctrl`
scores and labels:

| exon | score `s_i` | label `y_i` |
|---|---|---|
| enhanced_1 | `X[1, p]` | 1 |
| enhanced_2 | `X[2, p]` | 1 |
| … | … | … |
| control_1 | `X[k, p]` | 0 |
| control_2 | `X[k+1, p]` | 0 |
| … | … | … |

The AUC is computed via the **Mann–Whitney U identity** with mid-rank
tie handling:

```
AUC(p) = (R_pos(p) − n_pos * (n_pos + 1) / 2) / (n_pos * n_ctrl)
```

where `R_pos(p)` is the sum of mid-ranks of *positive-class* scores
in the pooled ranking at position `p`, `n_pos = n_cat`, and
`n_ctrl` is the control count. Equivalently:

```
AUC(p) = P(score(positive) > score(negative))
       + ½ * P(score(positive) = score(negative))
```

So `AUC(p) = 1` means every positive-class exon has a strictly
larger score at position `p` than every control exon; `AUC(p) = 0.5`
means the score is no better than random at that base; `AUC(p) = 0`
means perfect anti-correlation (regulated exons are *lower* than
control at that base).

Three properties that follow directly from this definition matter
for how you should read the line plot:

- **Tie handling.** Mid-ranks (`scipy.stats.rankdata(method='average')`)
  mean that when many exons have the same score (a *very* common
  case at a single base where most exons have zero crosslinks), the
  AUC sits near 0.5 rather than reading some arbitrary tie-breaking
  ordering as "discrimination". A position where everyone is zero
  gives AUC = 0.5 exactly.
- **Rank invariance.** Any strictly-monotone transform of the
  scores (z-scoring, log, exponentiation, etc.) leaves `AUC(p)`
  unchanged at *that single position*. The `--xl_score` choice
  matters because it changes the *relative ordering of scores
  across exons*, not because it changes the absolute values — see
  the "How `--xl_score` impacts these calculations" section.
- **Vectorised.** All positions are ranked at once
  (`scipy.stats.rankdata(matrix, axis=0)`), then the formula above
  is a single subtract-and-divide across the position axis. Cost
  is O(`n_exons * n_positions * log(n_exons)`); chr21 ~1000 exons
  × 601 positions runs in a few hundred milliseconds.

**On the plot.** The y-axis is the signed transform
`auc_signed(p) = 2 * (AUC(p) − 0.5)`, so 0 is the
no-information baseline and the sign matches the rest of the RNA
maps (positive = enriched, negative = depleted). The line is
optionally Gaussian-smoothed with `--smoothing` for readability;
the underlying `auc` column in the TSV is unsmoothed.

**Optional p-values.** With `--roc_n_perm > 0`, the
{category, control} labels are shuffled that many times. The
per-position p-value is the proportion of shuffles where
`|AUC_null(p) − 0.5| ≥ |AUC_obs(p) − 0.5|`, with a Laplace
correction of `(k+1)/(B+1)` to avoid `p = 0`. With `--roc_n_perm 0`
(default), the `pvalue` column is `NaN`.

### How the per-region ROC curve is calculated

The line plot answers "where does the score discriminate?" but a
canonical ROC curve answers "**how well does the score discriminate
overall?**". To get one number per exon, the per-exon row of the
coverage matrix is collapsed via the `--roc_aggregator` setting:

```
roc_aggregator = mean: s_i = (1 / n_pos) * sum_p X[i, p]
roc_aggregator = max:  s_i = max_p X[i, p]
roc_aggregator = both: compute both, write both to TSV, plot both
```

Then for each non-control category `c` vs. control, we sort all
`n_cat + n_ctrl` scores in decreasing order and sweep a threshold
through them, top to bottom:

```
At each threshold step k (the k-th-highest score):
    TP(k) = #{positive-class exons with score >= s_(k)}
    FP(k) = #{control exons       with score >= s_(k)}
    TPR(k) = TP(k) / n_pos
    FPR(k) = FP(k) / n_ctrl
```

Plotted as TPR vs FPR, that's the ROC curve. The curve always
passes through `(0, 0)` (no exons called positive) and `(1, 1)`
(all called positive); the line `y = x` is the no-information
baseline. Ties in scores collapse threshold steps so the curve
draws as a single segment rather than a staircase. AUC of the
curve is the same Mann–Whitney rank-AUC formula given above,
applied to the per-exon aggregate scores `s_i` (not per-position).

The legend prints the AUC for each category in each panel; you'll
typically see one panel per splice-site region (4 by default, 6
with `--all_sites`). Two PDFs are written when
`--roc_aggregator both`: `..._curves_mean.pdf` and
`..._curves_max.pdf`.

### When to read which output

- The **per-position AUC line plot** tells you *where* in the window
  the score is informative. A peak at, say, +50 nt downstream of the
  3'ss says "if you slice the score at this exact base, regulated
  exons rank above control there"; a flat 0 line says the score is
  uninformative at every base.
- The **per-region ROC curves** tell you *whether*, taken as a whole,
  the aggregated-over-the-window score is a useful predictor — and
  give you a single AUC number to report. The choice of aggregator
  matters: `mean` rewards broad coverage; `max` rewards sharp peaks.
  Compare them to learn whether your track is peaky or diffuse.

The two outputs are **not redundant**. Consider an AI track that
puts the same total probability mass in two different places: at
the 3'ss vs at -250 nt deep in the intron. The per-region ROC will
score them identically (same mean over the window); the
per-position AUC plot will trivially distinguish them.

### How `--xl_score` impacts these calculations

This is the lever for using BED column 5. The matrix entries `X[i, p]`
above are not literal counts — they're a function of `--xl_score`.

Below, "BED" means the `-x` file; "B.name" and "B.score" are columns 4
and 5 of the BED. "Within transcript `t`" means "among all overlap
rows where B.name = `t`". The rows are the matches produced by
`bedtools intersect -wa -wb` between the splice-site window BED and
the crosslink BED.

| `--xl_score` | What `X[i, p]` becomes | When to use |
|---|---|---|
| `ignore` (default) | Number of BED rows that overlap exon `i` at base `p`. Each match contributes `+1` regardless of B.score. | iCLIP / eCLIP integer count BEDs where each row already represents one observation. Bit-identical to legacy behaviour. |
| `raw` | Sum of `B.score` over BED rows overlapping exon `i` at base `p`. | The BED already carries the score you want to use, on a scale that's comparable across rows (e.g. eCLIP merged-crosslink counts, peak heights, pre-normalised AI scores). |
| `per_transcript_zscore` | Same as `raw`, but `B.score` is replaced by `(B.score − μ_t) / σ_t` *within each B.name = `t`* before the sum. `μ_t`, `σ_t` are the mean and ddof=0 std-dev of B.score within that transcript. Singletons (one entry per transcript) and zero-variance transcripts get a 0 score. | AI tracks whose absolute score magnitudes are **not comparable across transcripts** (e.g. each transcript has its own dynamic range; what counts as "high" in transcript A is different from transcript B). |
**Crucial subtlety for `per_transcript_zscore`.** The renormalisation
happens *within the overlap rows that actually intersect a splice-site
window*. If you slice your BED to only the entries near splice sites,
the within-transcript z-scoring operates on those slices only — which
is usually what you want. Pass the **full per-transcript BED** to `-x`
if you want the per-transcript baseline computed from the whole
transcript rather than from the slice.

**Which mode actually changes the ROC/AUC?** Because the per-position
AUC is rank-based, only changes that affect the **ordering of `s_i`
across exons** matter:

- `ignore` → `raw`: changes ordering whenever multiple BED rows
  overlap the same (exon, position) cell — usually a clear gain
  for eCLIP-like inputs where row scores differ. (This is the
  largest practical change you'll see; in the chr21 PTBP1 smoke
  test, silenced × middle_3ss × max AUC went from 0.44 → 0.66.)
- `raw` → `per_transcript_zscore`: changes ordering whenever the
  per-transcript mean differs between transcripts (i.e. high-signal
  transcripts get their bias removed). Reduces the effect of an
  exon being in a generally-high-signal transcript; sharpens
  position-specific contributions.
- Modes never change `ignore`'s row count of 0 (no overlap = 0
  score in every mode), so positions with no overlap give AUC =
  0.5 (all-tied) regardless of `--xl_score`.

If BED column 4 is `.` (placeholder) for every row,
`per_transcript_zscore` treats each overlap row as its own
singleton-transcript, which means **it collapses to `raw`**.
That's the right behaviour, but it also means a BED without a
meaningful `name` column gives you no benefit from the
per-transcript mode — make sure column 4 carries your transcript
ID before you reach for it.

The pipeline auto-detects non-trivial values in BED column 5 and
logs a `WARN` if you've left `--xl_score` at the default
`ignore` — you can pass it the more informative mode without
re-checking your BED.

**Sweeping modes in one run.** `--xl_score` accepts multiple values:

```
rnamaps -i ... -x ... -f ... -fi ... \
  --enrichment roc_auc \
  --xl_score ignore raw per_transcript_zscore \
  --roc_aggregator both --roc_n_perm 1000 -p AI
```

This computes the coverage matrices once per mode and writes
per-mode outputs (`{prefix}_RNAmap_roc_auc_xlscore-{mode}.{tsv,pdf}`),
while the xl_score-invariant outputs (heatmap, exon_length,
totalExonsCovered, categorised-exons TSV) are written once.

### Output columns for `roc_auc`

`{prefix}_RNAmap_roc_auc.tsv` — one row per (category, position, region):

| Column | Meaning |
|---|---|
| `name`, `position`, `label` | Category, 1-indexed position in the window, splice-site region. |
| `auc` | Raw AUC at this position (Mann–Whitney rank). `0.5` is no information. |
| `auc_signed` | `2 * (auc − 0.5)`. Sign matches enriched/depleted. |
| `auc_signed_smoothed` | `--smoothing`-Gaussian-smoothed version (plotted column). |
| `pvalue` | Per-position permutation p (NaN when `--roc_n_perm = 0`). |
| `coverage`, `number_exons`, `norm_coverage`, `control_coverage`, `control_number_exons`, `control_norm_coverage`, `fold_change` | Legacy parity columns. |

`{prefix}_RNAmap_roc_auc_region_auc.tsv` — headline summary, one row per (region, category, aggregator):

| Column | Meaning |
|---|---|
| `name`, `label`, `aggregator` | Category, region, `mean` or `max`. |
| `auc` | The single per-region ROC AUC for this combination. |
| `n_cat`, `n_ctrl` | Sample sizes used. |
| `pvalue` | Permutation p (NaN when `--roc_n_perm = 0`). |

`{prefix}_RNAmap_roc_auc_roc_curves.tsv` — FPR/TPR points for re-plotting, one row per threshold step.

---

## Shrinkage for log2fc

`log2fc` is great as a "how much more / less" summary, but it has a
well-known failure mode: when *both* coverage rates are near zero, a
single-exon flip can produce an enormous, unstable swing. A position
where the control has zero crosslinks and one of fifty enhanced exons
has a crosslink doesn't mean "infinite fold change" -- it means "we
don't have enough crosslinks here to say anything." Shrinkage is how
we ask the tool to be honest about that and pull those flaky values
toward zero.

`--shrinkage` selects how that's done. The default is `magnitude`,
which decides how much to shrink based purely on the **absolute size
of the coverage rates at each position**. It does *not* depend on
how many exons are in your category or control set.

### What each mode actually does

#### `magnitude` (default, recommended)

In one sentence: it asks "what is the *bigger* of the two rates at
this position? If it's well above 5%, that's real signal -- leave it
alone. If it's below 5%, both sides are sparse -- shrink the fold
change toward zero."

The intuition:

- Each position has a coverage **rate** = (number of exons with at
  least one crosslink at this base) / (number of exons in that
  category). It's just a fraction, like "37% of enhanced exons had
  a crosslink here."
- The "size" of the signal at this position is `max(cat_rate,
  ctrl_rate)`. This is the key choice: if either side has a
  substantial rate, we have real signal to talk about; if *both* are
  tiny, the fold change is unreliable no matter how big the ratio
  looks.
- Define a **shrinkage scale** `τ` (default `0.05`, i.e. 5% rate)
  and a **shrinkage weight**:

  ```
  weight = τ / (max(cat_rate, ctrl_rate) + τ)
  ```

  When `max` is much bigger than `τ`, the weight is close to 0 (no
  shrinkage). When `max` is much smaller than `τ`, the weight is
  close to 1 (heavy shrinkage).
- Then:

  ```
  shrunk_log2fc = (1 - weight) × log2((cat + ε) / (ctrl + ε))
  ```

  where `ε = 1e-3` is a tiny constant that keeps the log finite when
  one rate is exactly zero.

The result, in plain English:

- **Strong, well-supported peak (e.g. 80% cat, 10% ctrl):** `max = 0.8`,
  weight ≈ 0.06. The raw `log2(8) ≈ 3.0` becomes ~2.8. Preserved.
- **Real peak with zero control (e.g. 100% cat, 0% ctrl):** `max = 1.0`,
  weight ≈ 0.05. Still a clearly positive `log2fc ≈ 9` (instead of
  `+inf`). Preserved.
- **Single-exon flip (e.g. 2% cat, 0% ctrl):** `max = 0.02`, weight ≈
  0.71. Heavily shrunk: raw `log2(0.021/0.001) ≈ 4.4` becomes ~1.3.
- **Both rates tiny but ratio looks big (e.g. 0.5% cat, 0.1% ctrl):**
  `max = 0.005`, weight ≈ 0.91. Almost completely shrunk to zero.

The shrinkage strength is controlled by **one knob**:
`--shrinkage_scale τ`. The default `0.05` means "below 5% rate counts
as sparse." Use a smaller value (e.g. `0.01`) for gentler shrinkage
that only bites at very low rates, or a larger one (e.g. `0.1`) for
more aggressive shrinkage that touches moderate rates too.

#### `pseudocount` (legacy)

Adds a small flat number `ε` to both rates before the log:

```
log2fc(p) = log2((rate_cat(p) + ε) / (rate_ctrl(p) + ε))
```

By default `ε = max(1e-3, 0.01 × median_control_rate)` per region.
This is the previous default. It avoids the worst log-of-zero blow-ups
but the constant `ε` doesn't know whether the rates at a given
position are big or small -- so it under-tames the sparse cases.

Pass `--pseudocount 0.01` for a fixed value, or
`--pseudocount_frac 0.05` to make the adaptive epsilon larger.

#### `none`

No regularisation at all. Will produce `+inf` / `-inf` at positions
where one rate is exactly zero. Use only for diagnostics or when you
want to inspect the raw, unshrunken behaviour.

### Does the number of exons impact this metric?

**Short answer: no.** This was the explicit design goal. The
shrinkage decision is a pure function of the two rates at each
position. If you re-run on a subset of half your exons, the
`log2fc` point estimate doesn't move (only the bootstrap CI ribbon
widens, as it should -- you have less data, so you should be less
sure).

**Longer answer:**

- The raw `log2fc(p) = log2(rate_cat(p) / rate_ctrl(p))` is just a
  ratio of fractions. The *point estimate* is invariant to the number
  of exons: 30/100 vs 3/10 give the same fold change. What changes
  with `n` is *how much you should trust* that point estimate.
- The bootstrap confidence-interval ribbon already handles trust
  correctly. A fold change computed from a 50-exon category gives a
  wide ribbon; a fold change from 2000 exons gives a narrow ribbon.
  That's the right tool for "how sure am I about this number."
- The shrinkage's job is separate: it's about taming individual
  *point estimates* that look implausible because **both rates are
  tiny in absolute terms**. That's an absolute-rate question, not a
  sample-size question. So `magnitude` shrinkage looks only at the
  rate values.
- If you want to know whether a sparse-position log2fc is reliable,
  read off the CI ribbon. If the ribbon brackets zero, the bootstrap
  is already telling you the signal is consistent with noise --
  whether or not the shrinkage further pulled the central line down.

### The `delta` track is *not* shrunk

The `delta(p) = cat_rate(p) - ctrl_rate(p)` panel is the additive
contrast and is bounded in `[-1, 1]` by construction. It can't blow
up the way `log2fc` can, so no shrinkage is applied. Use `delta` for
the rate-scale view of the comparison; use `log2fc` for the
multiplicative, library-size-invariant view.

### When does shrinkage matter most?

- **Sparse coverage regions** (e.g. far from the splice site, or any
  region where most exons have no crosslink): big benefit. This is
  where raw `log2fc` is most misleading. The bigger of the two rates
  is small in absolute terms, so the magnitude weight kicks in.
- **Real peaks (high cat rate vs low ctrl, or vice versa):** no
  meaningful change -- the bigger rate is well above the shrinkage
  scale, so the contrast is preserved within a few percent of raw.
- **Both-sides-similar baseline regions:** equal rates give exactly
  `log2fc = 0` regardless of magnitude. Shrinkage cannot invent a
  difference where the rates agree.

### Output columns

When you use `bootstrap_contrast`, the output TSV
(`{prefix}_RNAmap_bootstrap_contrast.tsv`) carries these extra columns
so you can audit what shrinkage did:

| Column | Meaning |
|---|---|
| `shrinkage` | The mode used (`magnitude`, `pseudocount`, or `none`). |
| `tau` | The shrinkage scale used under `magnitude` mode (default 0.05). NaN under other modes. |
| `pseudocount` | The epsilon used under `--shrinkage pseudocount`. NaN under other modes. |

---

## Control-set hygiene

Even when categories were assigned by FDR, the "control" pool can
silently include genuinely regulated exons that fall below the detection
threshold. This biases every cat-vs-ctrl statistic toward null, regardless
of the test you use, and **bootstrapping does not fix this** — it's a
bias, not a variance. The `--control_set` flag offers two stricter
modes:

- `--control_set strict` keeps only control exons with
  `|dPSI| < --control_max_dpsi` (default 0.01) AND
  `FDR > --control_min_fdr` (default 0.5). Removes likely contaminated
  controls at the cost of a smaller (but cleaner) negative set.
- `--control_set constitutive_only` drops the original control pool and
  uses constitutive exons (high maxPSI, |dPSI| ≈ 0) as the negative set.
  By definition uninvolved in regulation, but typically fewer.

In VastDB mode, `strict` is partially applicable (no FDR column) and
will fall back to the dPSI cutoff with a warning. `constitutive_only`
works in both modes.

---

## Expression matching

When you provide `--gene_tpm`, RNA maps additionally filters the negative
set so that gene-expression levels are similar between regulated exons
(`enhanced` + `silenced`) and the comparison pool (`control` and, by
default, `constitutive`).

The current implementation uses quantile-stratified matching:

- It loads a gene-level TPM table (`gene_id`, `tpm`), auto-detects whether
  IDs overlap best as Ensembl IDs or symbols, and attaches TPM to each exon.
- It computes `log10(TPM + --tpm_pseudocount)` and bins regulated exons
  into `--tpm_n_bins` quantiles.
- It subsamples control (and constitutive unless
  `--no_match_constitutive`) to match those regulated-bin proportions as
  closely as possible.

This step runs after `--control_set` hygiene and before any enrichment
method, and is available in both rMATS and VastDB modes.

---

## Caveats

- The bootstrap and permutation tests do not control for covariates
  (exon length, GC content, expression). If you need covariate
  matching, pre-filter the relevant exons before running.
- For `bootstrap_contrast`, the central line and CI band reflect
  sampling uncertainty in the estimator, not contamination of the
  control pool. Use `--control_set` for the latter.
- For `cluster_perm`, the cluster p-value is bounded below by
  `1/(B+1)`, so very strong peaks all show `p ≈ 1/B`. Per-cluster
  effect size is still readable from the `mass` column.
- For `permutation_z`, the reported p comes from a normal
  approximation to the permutation null. The approximation is
  excellent in the body and very good in the tail for sums/means of
  dozens of exons, but extreme z values still carry Monte-Carlo noise
  in `μ_perm` and `σ_perm` — raise `--n_perm` if you need very stable
  scores at individual extreme positions.
- No multiple-testing correction is applied across positions outside of
  cluster_perm. Treat individual peaks as exploratory; biological
  replication of the curve shape is the strongest evidence.

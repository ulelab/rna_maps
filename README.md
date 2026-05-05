## RNA maps
Authors: charlotte.capitanchik@crick.ac.uk; leomwilkinson@gmail.com; aram.amalietti@gmail.com

---

## Overview

`rnamaps` generates RNA maps showing positional enrichment of RBP binding (from CLIP data) around regulated exons. It supports two input modes:

1. **rMATS mode** — takes rMATS differential splicing output and auto-categorises exons from dPSI/FDR thresholds
2. **VastDB mode** — takes pre-curated VastDB EVENT ID lists with categories already assigned

Both modes feed into the same analysis pipeline: splice site BED creation, CLIP coverage calculation, one or more enrichment analyses (bootstrap contrast / cluster permutation / legacy Fisher / legacy permutation z-score), RNA map plotting, per-exon heatmaps, and exon length distributions.

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
  -s, --smoothing       Smoothing window [DEFAULT: 15]
  --seed                Random seed for reproducible permutations / subsetting [DEFAULT: 42]
  -nc, --no_constitutive  Exclude constitutive category
  -ns, --no_subset      Disable subsetting of control/constitutive exons
                        (subsetting is auto-disabled for any non-fisher method)
  -ao, --all_sites      Include all 6 splice sites (default: 4 core sites)
  -p, --prefix          Prefix for output files
  --enrichment          One or more of {bootstrap_contrast, cluster_perm,
                        permutation_z, fisher} [DEFAULT: bootstrap_contrast].
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
                        negligible variance when n_ctrl >> n_cat.
  --pseudocount         Override adaptive log2FC pseudocount with a fixed
                        value [DEFAULT: adaptive].
  --pseudocount_frac    Adaptive pseudocount fraction of the regional
                        control coverage median [DEFAULT: 0.01].

Cluster-permutation options (--enrichment cluster_perm):
  --cluster_thresh      Cluster-defining |t| threshold [DEFAULT: 2.0]

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
| `bootstrap_contrast` *(default)* | per-position `delta` and `log2fc` with bootstrap CI ribbons | You want signed effect size + uncertainty, not a p-value. Especially good when n_ctrl >> n_cat or you want to compare maps from different libraries. |
| `cluster_perm` | per-cluster p-values, plot of Welch t with significance bars | You want significance per peak with FWER control across positions and don't want a saturating p-value floor. |
| `permutation_z` | signed -log10(p) and z-score per position | Legacy: per-position significance via label-permutation z-score. Useful for comparing to historical results. |
| `fisher` | signed -log10(p) per position | Legacy: per-position Fisher's exact test. Subsetting of control/constitutive is auto-enabled (others auto-disable it). |

You can pass multiple methods at once: `--enrichment bootstrap_contrast cluster_perm` produces the union of their outputs.

### `bootstrap_contrast` (default)

For each non-control category `c` vs. control:

1. Resample `n_c` exons with replacement from `c` and `n_ctrl` from
   control, `B = --n_boot` times.
2. Per iteration `b`, compute `mean_cov_c^b(p)` and `mean_cov_ctrl^b(p)`.
3. Report two contrasts per position with their `(2.5, 97.5)` percentile
   bands across `b`:

   - `delta(p) = mean_cov_c(p) - mean_cov_ctrl(p)` (additive scale).
   - `log2fc(p) = log2((mean_cov_c(p) + ε) / (mean_cov_ctrl(p) + ε))`
     (multiplicative scale, library-size invariant).

   Where the CI band excludes 0, the contrast is "significant" in a
   bootstrap sense, with no p-value involved.

**Class imbalance** (e.g. n_ctrl=10000, n_cat=500) is handled cleanly: the
bootstrap variance of `mean_cov_ctrl(p)` is negligible, so the CI on the
contrast reflects the cat-side variance, which is the right behaviour.
For a 5-10× speed-up at large `n_ctrl`, pass
`--bootstrap_control_fixed` to skip resampling control entirely (treat
its mean as a constant). The pipeline auto-suggests this when
`n_ctrl ≥ 20 × n_cat`.

**Pseudocount for log2fc** is adaptive by default:
`ε = max(1e-3, --pseudocount_frac × median(mean_cov_ctrl over region))`.
This avoids inflated log2FCs at sparse positions when `n_cat` is small.
Override with `--pseudocount FLOAT`.

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

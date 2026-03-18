## RNA maps
Authors: charlotte.capitanchik@crick.ac.uk; aram.amalietti@gmail.com; leomwilkinson@gmail.com

---

## Overview

`rna_maps.py` generates RNA maps showing positional enrichment of RBP binding (from CLIP data) around regulated exons. It supports two input modes:

1. **rMATS mode** — takes rMATS differential splicing output and auto-categorises exons from dPSI/FDR thresholds
2. **VastDB mode** — takes pre-curated VastDB EVENT ID lists with categories already assigned

Both modes feed into the same analysis pipeline: splice site BED creation, CLIP coverage calculation, Fisher's exact test enrichment, RNA map plotting, per-exon heatmaps, and exon length distributions.

---

## Quick Start

Create a conda environment and activate it:

```
conda env create -f environment.yml
conda activate rnamaps
```

### rMATS mode

```
python rna_maps.py \
  -i SE.MATS.JCEC.txt \
  -x CLIP_crosslinks.bed \
  -f genome.fa \
  -fi genome.fa.fai \
  -o output/ \
  -p PTBP1
```

### VastDB mode

```
python rna_maps.py \
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
python rna_maps.py -h

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
  --seed                Random seed for reproducible subsetting [DEFAULT: 42]
  -nc, --no_constitutive  Exclude constitutive category
  -ns, --no_subset      Disable subsetting of control/constitutive exons
  -ao, --all_sites      Include all 6 splice sites (default: 4 core sites)
  -p, --prefix          Prefix for output files

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
python rna_maps.py \
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
python=3.7.7
pandas=0.24.2
numpy=1.19.2
pybedtools=0.8.1
matplotlib=3.3.2
seaborn=0.11.0
scipy=1.3.1
```

---

## Reproducibility

The `--seed` flag (default: 42) controls the random seed used when subsetting control and constitutive exons to match the size of the largest regulated category. Setting the same seed produces identical subsets across runs. Disable subsetting entirely with `--no_subset`.
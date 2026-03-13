## RNA maps
Authors: charlotte.capitanchik@crick.ac.uk; aram.amalietti@gmail.com

---

> **Branch: VastDB_input**
> This branch introduces `rna_maps_vastdb_only.py`, a standalone script for generating RNA maps directly from VastDB EVENT IDs — no rMATS required. See [VastDB RNA Maps](#vastdb-rna-maps-rna_maps_vastdb_onlypy) below.

---

## VastDB RNA Maps (`rna_maps_vastdb_only.py`)

This script generates RNA maps from pre-categorised VastDB cassette exon EVENT IDs. Instead of running differential splicing analysis with rMATS, you supply lists of EVENT IDs already assigned to categories (enhanced, silenced, control, constitutive), along with a VastDB annotation file to look up genomic coordinates.

**Quick Start**

Create a conda environment and activate it:

```
conda env create -f environment.yml
conda activate rnamaps
```

Run with VastDB inputs:

```
python rna_maps_vastdb_only.py \
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

**Preparing inputs**

- **EVENT ID lists**: Plain text files with one VastDB EVENT ID per line (e.g. `HsaEX0012345`). Lines beginning with `#` are ignored.
- **VastDB annotation file**: The `EVENT_INFO-*.tab` file from VastDB, used to look up genomic coordinates (chromosome, exon start/end, flanking exon coordinates, strand) for each EVENT ID.
- **CLIP crosslinks**: iCLIP or eCLIP crosslinks in BED format.
- **Genome FASTA + index**: Reference genome `.fa` and `.fa.fai` files.

**Key differences from `rna_maps.py`**

| Feature | `rna_maps.py` | `rna_maps_vastdb_only.py` |
|---|---|---|
| Splicing input | rMATS `SE.MATS.JCEC.txt` | VastDB EVENT ID lists |
| Category assignment | Computed from dPSI / FDR thresholds | Pre-assigned by user |
| Annotation source | rMATS coordinates | VastDB `EVENT_INFO` file |
| rMATS dependency | Required | Not required |

**Usage**

```
python rna_maps_vastdb_only.py -h

required arguments:
  --vastdb_annotation   VastDB EVENT_INFO file
  -x, --clip            CLIP crosslinks BED file
  -f, --fasta           Genome FASTA file
  -fi, --fasta_index    Genome FASTA index (.fai)
  -o, --output_dir      Output directory
  -p, --prefix          Output file prefix

VastDB ID lists (at least one required):
  --vastdb_enhanced     Enhanced exon IDs (one per line)
  --vastdb_silenced     Silenced exon IDs (one per line)
  --vastdb_control      Control exon IDs (one per line)
  --vastdb_constitutive Constitutive exon IDs (one per line)

optional arguments:
  -w, --window          Window around splice sites [DEFAULT 300]
  -s, --smoothing       Smoothing window [DEFAULT 15]
  -nc, --no_constitutive  Exclude constitutive category
  -ns, --no_subset      Disable subsetting of control/constitutive to match enhanced/silenced counts
  -ao, --all_sites      Include all 6 splice sites (default: 4 core sites)
```

**Outputs**

- `{prefix}_VastDB_with_categories.tsv` — categorised exons with coordinates
- `{prefix}_RNAmap.tsv` — per-position coverage and enrichment data
- `{prefix}_RNAmap_-log10pvalue.pdf` — RNA map plot
- `execution_*.log` — run log

---

## Original RNA Maps (`rna_maps.py`)

The original script accepts rMATS quantified splicing files and assigns categories based on dPSI and FDR thresholds.

**Quick Start**

```
python rna_maps.py \
-i test/chr21_PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt \
-x test/chr21_hela_ptbp1_iclip_sorted_merged.bed \
-f test/homosapien-hg37-chr21.fa \
-fi test/homosapien-hg37-chr21.fa.fai
```

**Preparing RNA-Seq data**:

This code accepts rMATs quantified files for cassette exons (e.g. SE.MATS.JCEC.txt).

If your condition is RBP knockdown be sure to run your comparison as condition - control, such that definitions of enhanced and repressed are correct. If your condition is RBP overexpression you will need to run the comparison as control - condition. In the generic example group1 - group2 consider that the definition of "enhanced" or "repressed" are in reference to group2. ie. an exon is enhanced in group2 vs. group1.

**Multivalency analysis**:

Multivalency analysis adds on run time & involves installing the Ule lab's GeRMs package which is still in development, so it is optional and enabled with the flag `-v`.
Currently to run the analysis you will need to install the GeRMs package. To do this clone the repository to your computer somewhere and run the following command from within the repository (you will need to have R devtools installed):
`R -e 'devtools::install()'`

You will need to ensure you have the GeRMs requirements installed too, which are: biostrings, parallel, logger and optparse.
Finally, when you run RNA maps you will need to provide the location of your "germs" repo, so that the script can find germs.R to run the multivalency calculations using the flag `-g`, so our test command for running multivalency will look like:
```
python rna_maps.py \
-i test/chr21_PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt \
-x test/chr21_hela_ptbp1_iclip_sorted_merged.bed \
-f test/homosapien-hg37-chr21.fa \
-fi test/homosapien-hg37-chr21.fa.fai \
-v -g ../germs
```
If you want to create a multivalency map alone (no CLIP data) simply run the above command with the `-x crosslinks.bed` excluded.

**Usage**:
```
python rna_maps.py -h
usage: rna_maps.py [-h] -i INPUTSPLICE [-x [INPUTXLSITES]] -f GENOMEFASTA -fi
                   FASTAINDEX [-o [OUTPUTPATH]] [-w [WINDOW]] [-s [SMOOTHING]]
                   [-mc [MINCTRL]] [-xc [MAXCTRL]] [-xi [MAXINCL]] [-xf [MAXFDR]]
                   [-xe [MAXENH]] [-ms [MINSIL]] [-v] [-nc] [-ns] [-ao] [-g [GERMSDIR]]
                   [-p PREFIX]

Plot CLIP crosslinks around regulated exons to study position-dependent impact on pre-
mRNA splicing.

required arguments:
  -i INPUTSPLICE, --inputsplice INPUTSPLICE
                        quantification of differential splicing produced by rMATS
  -f GENOMEFASTA, --genomefasta GENOMEFASTA
                        genome fasta file (.fa)
  -fi FASTAINDEX, --fastaindex FASTAINDEX
                        genome fasta index file (.fai)

options:
  -h, --help            show this help message and exit
  -x [INPUTXLSITES], --inputxlsites [INPUTXLSITES]
                        CLIP crosslinks in BED file format
  -o [OUTPUTPATH], --outputpath [OUTPUTPATH]
                        output folder [DEFAULT current directory]
  -w [WINDOW], --window [WINDOW]
                        window around regulated splicing events to plot crosslinks
                        [DEFAULT 300]
  -s [SMOOTHING], --smoothing [SMOOTHING]
                        smoothing window for plotting crosslink signal [DEFAULT 15]
  -mc [MINCTRL], --minctrl [MINCTRL]
                        minimum dPSI for control events [DEFAULT -0.05]
  -xc [MAXCTRL], --maxctrl [MAXCTRL]
                        maximum dPSI for control events [DEFAULT 0.05]
  -xi [MAXINCL], --maxincl [MAXINCL]
                        maximum PSI for control exons, above this limit exons are
                        considered constitutive [DEFAULT 0.9]
  -xf [MAXFDR], --maxfdr [MAXFDR]
                        maximum FDR for regulated events, above this events fall in
                        "rest" class, is used for rMATS [DEFAULT 0.1]
  -xe [MAXENH], --maxenh [MAXENH]
                        maximum inclusion for exons to be considered enhanced [DEFAULT
                        -0.05]
  -ms [MINSIL], --minsil [MINSIL]
                        minimum inclusion for exons to be considered silenced [DEFAULT
                        0.05]
  -v, --multivalency
  -nc, --no_constitutive
                        Exclude constitutive category from the output
  -ns, --no_subset      Disable subsetting of control/constitutive exons to match
                        enhanced/silenced counts
  -ao, --all_sites      Include all splice sites (upstream_3ss and downstream_5ss),
                        default is core sites only
  -g [GERMSDIR], --germsdir [GERMSDIR]
                        directory for where to find germs.R for multivalency analysis
                        eg. /Users/Bellinda/repos/germs [DEFAULT current directory]
  -p PREFIX, --prefix PREFIX
                        prefix for output files [DEFAULT inputsplice file name]
```

**Dependencies** (these are the versions the script was developped with, pandas >= 1 introduced breaking changes, please use these versions):
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

### Definitions

**Event types**

Control: An event that doesn't change in inclusion (PSI) in this RBP knockdown, but might in another circumstance. Typical definition:

```
dPSI   ( -1 <---------- - 0.05xxxxx0xxxxx0.05----------> 1 )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9------> 1 )
FDR    (  0 xxxxx0.1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 1 )
```

Constituitive: An event that doesn't change in inclusion (PSI) in this RBP knockdown, but is unlikely to change in another circumstance. Typically defined as a control event *plus* have a maximum inclusion (PSI) of > 0.9-0.99.

```
dPSI   ( -1 <---------- - 0.05xxxxx0xxxxx0.05----------> 1 )
maxPSI (  0 ----------------------------------0.9xxxxxx> 1 )
FDR    (  0 xxxxx0.1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 1 )
```

Enhanced: An event that is *less* included in RBP knockdown, suggesting the RBP *promotes/enhances* inclusion of the event.

```
dPSI   ( -1 <xxxxxxxxxxx- 0.05-----0-----0.05----------> 1 )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9xxxxxx> 1 )
FDR    (  0 xxxxx0.1-----------------------------------> 1 )
```

Silenced: An event that is *more* included in RBP knockdown, suggesting the RBP *represses/silences* inclusion of the event.

```
dPSI   ( -1 <---------- - 0.05-----0-----0.05xxxxxxxxxx> 1 )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9xxxxxx> 1 )
FDR    (  0 xxxxx0.1-----------------------------------> 1 )
```

Enhanced/Silenced rest: A silenced or enhanced event where the FDR does not fall below the threshold.

```
dPSI   (            As in silenced or enhanced             )
maxPSI (  0 xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx0.9xxxxxx> 1 )
FDR    (  0 -----0.1xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx> 1 )
```

**Hierarchy**

When it comes to alternative exons, an exon may be involved in multiple events, but we want to avoid plotting it many times, so we implement a hierarchy:

1. If an exon meets criteria for silenced or enhanced this is designated, if criteria for both is met the most extreme dPSI value is preferred.
2. Of remaining exons, if they meet criteria for enhanced/silenced rest this is designated, if criteria for both is met the most extreme dPSI value is preferred.
3. Of remaining exons, if they meet critera for constituitive, this is designated.
4. Of remaining exons, if they meet critera for control, this is designated.

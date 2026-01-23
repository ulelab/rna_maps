## RNA maps
Authors: charlotte.capitanchik@crick.ac.uk; aram.amalietti@gmail.com

**Quick Start**

Create a conda environment with all dependencies and activate it:

```
conda env create -f environment.yml
conda activate rnamaps
```

Then run the test dataset to check the code is working:

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

**Usage**:  
```
python rna_maps.py -h
usage: rna_maps.py [-h] -i INPUTSPLICE [-x [INPUTXLSITES]] -f GENOMEFASTA -fi
                   FASTAINDEX [-o [OUTPUTPATH]] [-w [WINDOW]] [-s [SMOOTHING]]
                   [-mc [MINCTRL]] [-xc [MAXCTRL]] [-xi [MAXINCL]] [-xf [MAXFDR]]
                   [-xe [MAXENH]] [-ms [MINSIL]] [-v] [-nc] [-ns] [-g [GERMSDIR]]
                   [-p PREFIX]

Plot CLIP crosslinks around regulated exons to study position-dependent impact on pre-mRNA splicing.

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
  -g [GERMSDIR], --germsdir [GERMSDIR]
                        directory for where to find germs.R for multivalency analysis
                        eg. /Users/Bellinda/repos/germs [DEFAULT current directory]
  -p PREFIX, --prefix PREFIX
                        prefix for output files [DEFAULT inputsplice file name
```

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

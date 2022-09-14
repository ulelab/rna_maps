## RNA maps
Author: aram.amalietti@gmail.com

**Quick Start**

Create a conda environment with all dependencies and activate it:

```
conda env create -f environment.yml
conda activate rnamaps
```

Then run the test dataset to check the code is working:

```
python3 \
rna_maps.py \
test/chr21_PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt \
test/chr21_hela_ptbp1_iclip_sorted_merged.bed \
test/GRCh38.release34.primary_assembly.genome.fa.fai
```
*Preparing RNA-Seq data*:

This code accepts rMATs junction only quantified files, or Whippet.
Be sure to run your comparison as condition - control, such that definitions of enhanced and repressed are correct.


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
python3 rna_maps.py -h                                                                                                                                                      
usage: rna_maps.py [-h] -i INPUTSPLICE -x INPUTXLSITES -f FASTAINDEX
                   [-o [OUTPUTPATH]] [-w [WINDOW]] [-s [SMOOTHING]]
                   [-mc [MINCTRL]] [-xc [MAXCTRL]] [-xi [MAXINCL]]
                   [-xf [MAXFDR]] [-xe [MAXENH]] [-ms [MINSIL]]

Plot CLIP crosslinks around regulated exons to study position-dependent impact
on pre-mRNA splicing.

required arguments:
  -i INPUTSPLICE, --inputsplice INPUTSPLICE
                        quantification of differential splicing produced by
                        rMATS
  -x INPUTXLSITES, --inputxlsites INPUTXLSITES
                        CLIP crosslinks in BED file format
  -f FASTAINDEX, --fastaindex FASTAINDEX
                        genome fasta index file (.fai)

optional arguments:
  -h, --help            show this help message and exit
  -o [OUTPUTPATH], --outputpath [OUTPUTPATH]
                        output folder [DEFAULT current directory]
  -w [WINDOW], --window [WINDOW]
                        window around regulated splicing events to plot
                        crosslinks [DEFAULT 300]
  -s [SMOOTHING], --smoothing [SMOOTHING]
                        smoothing window for plotting crosslink signal
                        [DEFAULT 15]
  -mc [MINCTRL], --minctrl [MINCTRL]
                        minimum dPSI for control events [DEFAULT -0.05]
  -xc [MAXCTRL], --maxctrl [MAXCTRL]
                        maximum dPSI for control events [DEFAULT 0.05]
  -xi [MAXINCL], --maxincl [MAXINCL]
                        maximum PSI for control exons, above this limit exons
                        are considered constitutive [DEFAULT 0.9]
  -xf [MAXFDR], --maxfdr [MAXFDR]
                        maximum FDR for regulated events, above this events
                        fall in "rest" class, is used for rMATS [DEFAULT 0.1]
  -xe [MAXENH], --maxenh [MAXENH]
                        maximum inclusion for exons to be considered enhanced
                        [DEFAULT -0.05]
  -ms [MINSIL], --minsil [MINSIL]
                        minimum inclusion for exons to be considered silenced
                        [DEFAULT 0.05] 
```

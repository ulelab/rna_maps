#!/bin/bash

sbatch --mem 64GB --time=0:30:00 --wrap="\
python3 \
../rna_maps.py \
chr21_PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt \
chr21_hela_ptbp1_iclip_sorted_merged.bed \
GRCh38.release34.primary_assembly.genome.fa.fai \
. \
1000 \
15 \
-0.05 \
0.05 \
0.9 \
0.1 \
-0.05 \
0.05 \
0.9"

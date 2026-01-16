#!/bin/bash

python3 ~/repos/rna_maps/rna_maps.py \
-i ~/repos/rna_maps/PTBP1_2_Gueroussov2015_SE.MATS.JCEC.txt \
-x PTBP1_iCLIP.xl.clusters.bed \
-f ~/data/ref/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
-fi ~/data/ref/GRCh38.primary_assembly.genome.fa.fai \
-p "PTBP1_iCLIP.xl.clusters"
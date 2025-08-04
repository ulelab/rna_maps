#!/bin/bash

python3 ../rna_maps.py \
-i ../RBFOX2_encodeHEPG2_SE.MATS.JCEC.txt \
-x ../HepG2_RBFOX2_clippy_rollmean50_stdev1.0_minGeneCount5_broadPeaks.bed \
-f ../../../data/ref/GRCh38.primary_assembly.genome.fa \
-fi ../../../data/ref/GRCh38.primary_assembly.genome.fa.fai \
-w 300
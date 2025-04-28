#!/bin/bash

python3 ../rna_maps.py \
-i ../Karen_12hr_PPIG_SE.MATS.JCEC.txt \
-f ../../../data/ref/GRCh38.primary_assembly.genome.fa \
-fi ../../../data/ref/GRCh38.primary_assembly.genome.fa.fai \
-v -g ../../germ \
-xf 0.01 -ms 0.4 -xe -0.4
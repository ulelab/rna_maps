## rna_maps
Author: aram.amalietti@gmail.com


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
**Preparing RNA-Seq data**:

This code accepts rMATs junction only quantified files, or Whippet.
Be sure to run your comparison as condition - control, such that definitions of enhanced and repressed are correct.

**Usage**:  
```
    python3 <path_to_script> <de_file_path> <sites_file_path> <genome_fai_path> <output_folder> <window> <smoothing> <min_ctrl> <max_ctrl> <max_inclusion> <max_fdr> <max_enc> <min_sil> <min_prob_whippet>
```
`de_file`:*file with differential splicing table, rMATS and whippet are supported;*  
`xl_bed`:*BED file with genomic coordinates of landmarks that are used for mapping around exons;*  
`fai`:*FASTA index file;*  
`output_folder`:*folder where the results will be saved, make sure it exists and is writable;*  
`window`:*flanks around exons splice sites where coverages are mapped(recommended 300);*  
`smoothing`:*size of smoothing window(recommended 15);*   
`min_ctrl`:*minimal inclusion change for control exons(recommended -0.05);*  
`max_ctrl`:*maximal inclusion change for control exons(recommended 0.05);*  
`max_inclusion`:*maximal inclusion for control exons, above this limit exons are considered constitutive(recommended 0.9 to 0.99);*  
`max_fdr`:*maximal FDR for regulated exons, above exons fall in rest class, is used for rMATS(recommended 0.1);*  
`max_enc`:*maximum inclusion for exons to be considered enhanced (recommended -0.05);*  
`min_sil`:*minimum inclusion for exons to be considered silenced (recommended 0.05);*  
`min_prob_whippet`:*minimum probability for exons to considered regulated, below exons are labeled as rest, is used for whippet (recommended 0.9);*  

#!/bin/bash

# Pipeline for analysing environmental samples.
# Part 4 - Making analyses
# Requires rep_set_final.fna, otu_table_final.txt and a mapping file.

echo "Launching....."

# Only need if running macqiime
source /macqiime/configs/bash_profile.txt

# Make_Tree
mkdir 11_Tree
mkdir 11_Tree/16S
mkdir 11_Tree/18S
mkdir 11_Tree/UPA
mkdir 11_Tree/tufA

make_phylogeny.py -i 09_Alignment/16S/rep_set_aligned_biomfiltered.fasta -o 11_Tree/16S/tree.nwk
make_phylogeny.py -i 09_Alignment/18S/rep_set_aligned_biomfiltered.fasta -o 11_Tree/18S/tree.nwk
make_phylogeny.py -i 09_Alignment/UPA/rep_set_aligned_biomfiltered.fasta -o 11_Tree/UPA/tree.nwk
make_phylogeny.py -i 09_Alignment/tufA/rep_set_aligned_biomfiltered.fasta -o 11_Tree/tufA/tree.nwk

# Core diversity analyses - depth = 10000 (100 for tufA while alignment isn't working)
mkdir 12_Core_Diversity
core_diversity_analyses.py -i 10_OTU_filtering/16S/otu_table_final.biom -m 05_Mapping/16S/Map_16S.txt -t 11_Tree/16S/tree.nwk -e 10000 -o 12_Core_Diversity/16S/
core_diversity_analyses.py -i 10_OTU_filtering/18S/otu_table_final.biom -m 05_Mapping/18S/Map_18S.txt -t 11_Tree/18S/tree.nwk -e 10000 -c "Patchiness_host,Patchiness_sample" -o 12_Core_Diversity/18S/
core_diversity_analyses.py -i 10_OTU_filtering/UPA/otu_table_final.biom -m 05_Mapping/UPA/Map_UPA.txt -t 11_Tree/UPA/tree.nwk -e 10000 -c "Patchiness_host,Patchiness_sample" -o 12_Core_Diversity/UPA/
core_diversity_analyses.py -i 10_OTU_filtering/tufA/otu_table_final.biom -m 05_Mapping/tufA/Map_tufA.txt -t 11_Tree/tufA/tree.nwk -e 10000 -c "Patchiness_host,Patchiness_sample" -o 12_Core_Diversity/tufA/


#!/bin/bash

# Pipeline for analysing environmental samples.
# Beta diversity analysis
# Requires rep_set_final.fna, otu_table_final.txt and a mapping file.

echo "Launching....."


# Make_Tree
echo "Making phylogeny"
mkdir 11_Tree
make_phylogeny.py -i 09_Euk_Otus/16S.euk_otus_aln.fas -o 11_Tree/16S_tree.nwk
make_phylogeny.py -i 09_Euk_Otus/18S.euk_otus_aln.fas -o 11_Tree/18S_tree.nwk
make_phylogeny.py -i 09_Euk_Otus/UPA.euk_otus_aln.fas -o 11_Tree/UPA_tree.nwk
make_phylogeny.py -i 09_Euk_Otus/tufA.euk_otus_aln.fas -o 11_Tree/tufA_tree.nwk


# Core diversity analyses - depth = 200.
# Here - Eukaryotic OTUs only. For beta diversity, better use all seqs.
echo "core diversity analysis of Eukariotic OTUs - check seq. depth!"
mkdir 12_Core_Diversity
core_diversity_analyses.py -i 07_OTU_table/16S/OTU_table_6.biom -m Mapping/map.txt -t 11_Tree/16S_tree.nwk -e 1000 -o 12_Core_Diversity/16S/
core_diversity_analyses.py -i 07_OTU_table/18S/OTU_table_6.biom -m Mapping/map.txt -t 11_Tree/18S_tree.nwk -e 1000 -o 12_Core_Diversity/18S/
core_diversity_analyses.py -i 07_OTU_table/UPA/OTU_table_6.biom -m Mapping/map.txt -t 11_Tree/UPA_tree.nwk -e 1000 -o 12_Core_Diversity/UPA/
core_diversity_analyses.py -i 07_OTU_table/tufA/OTU_table_6.biom -m Mapping/map.txt -t 11_Tree/tufA_tree.nwk -e 1000 -o 12_Core_Diversity/tufA/

echo "Done"
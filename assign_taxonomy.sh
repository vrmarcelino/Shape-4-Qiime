#!/bin/bash

#PBS -V
#PBS -N Pipe_tufA
#PBS -k eo
#PBS -l nodes=3:ppn=3
# I think the PBS stuff isn't working...

### Assign taxonomy to the OTUs.
### Using here RDP. You may prefer to use uclust for 16S and 18S, but it doesn't work well for UPA and tufA.
### Input files: taxonomy map and fasta file of reference rep_sequences
### OTUs.fasta from blast filter



echo "Assign taxomony - RDP - Qiime"
mkdir 06_Taxonomy

assign_taxonomy.py -i 04_blastfilter/02_final_results/16S.otus.fasta -t ~/EnvSeq_reference_datasets/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -r ~/EnvSeq_reference_datasets/gg_13_8_otus/rep_set/97_otus.fasta -o 06_Taxonomy/16S -m rdp >> screen_16S.out 2>> screen_16S.err
assign_taxonomy.py -i 04_blastfilter/02_final_results/18S.otus.fasta -t ~/EnvSeq_reference_datasets/Silva_111/taxonomy/97_Silva_111_taxa_map_RDP_6_levels.txt -r ~/EnvSeq_reference_datasets/Silva_111/rep_set/97_Silva_111_rep_set.fasta -o 06_Taxonomy/18S -m rdp >> screen_18S.out 2>> screen_18S.err
assign_taxonomy.py -i 04_blastfilter/02_final_results/UPA.otus.fasta -t ~/EnvSeq_reference_datasets/Ref_dataset_UPA/id_to_taxonomy_map.txt -r ~/EnvSeq_reference_datasets/Ref_dataset_UPA/reference_sequences.fna -o 06_Taxonomy/UPA -m rdp >> screen_UPA.out 2>> screen_UPA.err
assign_taxonomy.py -i 04_blastfilter/02_final_results/tufA.otus.fasta -t ~/EnvSeq_reference_datasets/Ref_dataset_tufA/id_to_taxonomy_map.txt -r ~/EnvSeq_reference_datasets/Ref_dataset_tufA/reference_sequences.fna -o 06_Taxonomy/tufA -m rdp >> screen_tufA.out 2>> screen_tufA.err

echo "done"




#!/bin/bash
#PBS -V
#PBS -N Pipe2_QIIME
#PBS -k eo
#PBS -l nodes=3:ppn=3

# Pipeline for analysing environmental samples.
# Part 2 - Quality control and picking OTUs
# Works with Qiime 1.9
# 19 - Jan - 2015

# Note: Seems that Qiime 1.9 doesn't make html OTU tables and html heatmaps anymore. Why???
# Use the tax assignments and do it in Macqiime v. 1.8 to check contamination rates

echo "Launching....."


### Quality control using Split Library
# Quality filter: min length (-l) 150bp for 16S and 200 for other genes, quality score (-s = 25) default = 25, max-homopolymer = 10 for UPA, 6 (default) for others
echo "Split_libraries running"
split_libraries.py -m 05_Mapping/16S/Map_16S.txt -f 04_All_records/16S/all_records.fna -q 04_All_records/16S/all_records.qual -p -s 25 -l 150 -b 16 -o 06_Split_Library/16S >> screen_16S.out 2>> screen_16S.err
split_libraries.py -m 05_Mapping/18S/Map_18S.txt -f 04_All_records/18S/all_records.fna -q 04_All_records/18S/all_records.qual -p -s 25 -l 200 -b 16 -o 06_Split_Library/18S >> screen_18S.out 2>> screen_18S.err
split_libraries.py -m 05_Mapping/UPA/Map_UPA.txt -f 04_All_records/UPA/all_records.fna -q 04_All_records/UPA/all_records.qual -p -s 25 -l 200 -b 16 -H 10 -o 06_Split_Library/UPA >> screen_UPA.out 2>> screen_UPA.err
split_libraries.py -m 05_Mapping/tufA/Map_tufA.txt -f 04_All_records/tufA/all_records.fna -q 04_All_records/tufA/all_records.qual -p -s 25 -l 200 -b 16 -o 06_Split_Library/tufA >> screen_tufA.out 2>> screen_tufA.err


### Chimera checking
echo "Chimera checking with usearch 6.1 ..."
mkdir 07_Chimeras
identify_chimeric_seqs.py -i 06_Split_Library/16S/seqs.fna -m usearch61 -o 07_Chimeras/16S/ -r gg_13_8_otus/rep_set/97_otus.fasta >> screen_16S.out 2>> screen_16S.err
filter_fasta.py -f 06_Split_Library/16S/seqs.fna -o 07_Chimeras/16S/seqs_no_chimeras.fna -s 07_Chimeras/16S/chimeras.txt -n >> screen_16S.out 2>> screen_16S.err

identify_chimeric_seqs.py -i 06_Split_Library/18S/seqs.fna -m usearch61 -o 07_Chimeras/18S/ -r Silva_111/rep_set/97_Silva_111_rep_set.fasta >> screen_18S.out 2>> screen_18S.err
filter_fasta.py -f 06_Split_Library/18S/seqs.fna -o 07_Chimeras/18S/seqs_no_chimeras.fna -s 07_Chimeras/18S/chimeras.txt -n >> screen_18S.out 2>> screen_18S.err

identify_chimeric_seqs.py -i 06_Split_Library/UPA/seqs.fna -m usearch61 -o 07_Chimeras/UPA/ -r Ref_dataset_UPA/reference_sequences.fna >> screen_UPA.out 2>> screen_UPA.err
filter_fasta.py -f 06_Split_Library/UPA/seqs.fna -o 07_Chimeras/UPA/seqs_no_chimeras.fna -s 07_Chimeras/UPA/chimeras.txt -n >> screen_UPA.out 2>> screen_UPA.err

identify_chimeric_seqs.py -i 06_Split_Library/tufA/seqs.fna -m usearch61 -o 07_Chimeras/tufA/ -r Ref_dataset_tufA/reference_sequences.fna >> screen_tufA.out 2>> screen_tufA.err
filter_fasta.py -f 06_Split_Library/tufA/seqs.fna -o 07_Chimeras/tufA/seqs_no_chimeras.fna -s 07_Chimeras/tufA/chimeras.txt -n >> screen_tufA.out 2>> screen_tufA.err


### Pick OTUs, make OTU table and align
echo "Picking OTUs...it takes a while..."
# -a and -O for running in parallel (3 cores), --min_otu_size (to retain OTU) = 5
mkdir 09_Alignment

#16S follows the qiime pipeline - assigns tax. with Uclust
echo "Picking open reference 16S OTUs"
pick_open_reference_otus.py -i 07_Chimeras/16S/seqs_no_chimeras.fna -o 08_OTUs/16S -r gg_13_8_otus/rep_set/97_otus.fasta -a -O 4 --min_otu_size 5 >> screen_16S.out 2>> screen_16S.err
mkdir 09_Alignment/16S
mv 08_OTUs/16S/pynast_aligned_seqs/* 09_Alignment/16S
rm -r 08_OTUs/16S/pynast_aligned_seqs/

# 18S Can't align or assign taxonomy, so these will be done step by step - assigns tax. with Uclust
echo "Picking open reference 18S OTUs" 
pick_open_reference_otus.py -i 07_Chimeras/18S/seqs_no_chimeras.fna -o 08_OTUs/18S -r Silva_111/rep_set/97_Silva_111_rep_set.fasta --suppress_taxonomy_assignment --suppress_align_and_tree -a -O 4 --min_otu_size 5 >> screen_18S.out 2>> screen_18S.err
assign_taxonomy.py -i 08_OTUs/18S/rep_set.fna -t Silva_111/taxonomy/97_Silva_111_taxa_map.txt -r Silva_111/rep_set/97_Silva_111_rep_set.fasta -o 08_OTUs/18S/uclust_assigned_taxonomy >> screen_18S.out 2>> screen_18S.err
align_seqs.py -i 08_OTUs/18S/rep_set.fna -t Silva_111/rep_set_aligned/97_Silva_111_rep_set.fasta -o 09_Alignment/18S/ >> screen_18S.out 2>> screen_18S.err
filter_alignment.py -i 09_Alignment/18S/rep_set_aligned.fasta -s -o 09_Alignment/18S/ >> screen_18S.out 2>> screen_18S.err


# UPA - de novo OTU picking. Assigns tax. with RDP classifier. Can't use --min_count to remove OTUs with less than 5 reads here. Take them out later.
echo "Picking de novo UPA OTUs"
#0.90 similarity treshould
pick_otus.py -i 07_Chimeras/UPA/seqs_no_chimeras.fna -s 0.90 -o 08_OTUs/UPA --threads 4 >> screen_UPA.out 2>> screen_UPA.err
pick_rep_set.py -i 08_OTUs/UPA/seqs_no_chimeras_otus.txt -f 07_Chimeras/UPA/seqs_no_chimeras.fna -m most_abundant -o 08_OTUs/UPA/rep_set.fna >> screen_UPA.out 2>> screen_UPA.err
assign_taxonomy.py -i 08_OTUs/UPA/rep_set.fna -t Ref_dataset_UPA/id_to_taxonomy_map.txt -r Ref_dataset_UPA/reference_sequences.fna -o 08_OTUs/UPA/rdp_assigned_taxonomy -m rdp >> screen_UPA.out 2>> screen_UPA.err
align_seqs.py -i 08_OTUs/UPA/rep_set.fna -t Ref_dataset_UPA/reference_sequences_aligned.fna -o 09_Alignment/UPA/ >> screen_UPA.out 2>> screen_UPA.err
filter_alignment.py -i 09_Alignment/UPA/rep_set_aligned.fasta -s -o 09_Alignment/UPA/ >> screen_UPA.out 2>> screen_UPA.err


# tufA - de novo OTU picking. Assigns tax. with RDP classifier. Can't use --min_count to remove OTUs with less than 5 reads here. Take them out later.
echo "Picking de novo tufA OTUs"
#0.85 similarity treshould
pick_otus.py -i 07_Chimeras/tufA/seqs_no_chimeras.fna -s 0.85 -o 08_OTUs/tufA --threads 4 >> screen_tufA.out 2>> screen_tufA.err
pick_rep_set.py -i 08_OTUs/tufA/seqs_no_chimeras_otus.txt -f 07_Chimeras/tufA/seqs_no_chimeras.fna -m most_abundant -o 08_OTUs/tufA/rep_set.fna >> screen_tufA.out 2>> screen_tufA.err
assign_taxonomy.py -i 08_OTUs/tufA/rep_set.fna -t Ref_dataset_tufA/id_to_taxonomy_map.txt -r Ref_dataset_tufA/reference_sequences.fna -o 08_OTUs/tufA/rdp_assigned_taxonomy -m rdp >> screen_tufA.out 2>> screen_tufA.err
align_seqs.py -i 08_OTUs/tufA/rep_set.fna -t Ref_dataset_tufA/reference_sequences_aligned.fna -o 09_Alignment/tufA/ >> screen_tufA.out 2>> screen_tufA.err
filter_alignment.py -i 09_Alignment/tufA/rep_set_aligned.fasta -s -o 09_Alignment/tufA/ >> screen_tufA.out 2>> screen_tufA.err


echo "End of pipe 2"
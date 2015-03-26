#!/bin/bash

#PBS -V
#PBS -N Pipe_tufA
#PBS -k eo
#PBS -l nodes=3:ppn=3

# Analyse tufA of all runs (1,2,3) merged

# Run separate_genes_w_global_trimming.py before!! 

# At species level (98%) get a number of OTUs and get a tree.

### Split library using the -n option to rename the seqs according run. Increased minimum lenght to 400 (ensure the global alignment thing)
# Split libraries here is mainly to add labels to the sequences, the quality control will be done using UPARSE.
# parameters: min length (-l) 400bp , quality score (-s = 1), max-homopolymer = 10. Save qual file.
echo "Split_libraries running - Qiime"
split_libraries.py -m 05_Mapping/R1_Map_tufA.txt -f 04_All_records/R1_all_records.fna -q 04_All_records/R1_all_records.qual -p -s 1 -l 400 -b 16 -H 10 -n 10000000 --record_qual_scores -o 06_Split_Library/run1 >> screen_tufA.out 2>> screen_tufA.err
split_libraries.py -m 05_Mapping/R2_Map_tufA.txt -f 04_All_records/R2_all_records.fna -q 04_All_records/R2_all_records.qual -p -s 1 -l 400 -b 16 -H 10 -n 20000000 --record_qual_scores -o 06_Split_Library/run2 >> screen_tufA.out 2>> screen_tufA.err
split_libraries.py -m 05_Mapping/R3_Map_tufA.txt -f 04_All_records/R3_all_records.fna -q 04_All_records/R3_all_records.qual -p -s 1 -l 400 -b 16 -H 10 -n 30000000 --record_qual_scores -o 06_Split_Library/run3 >> screen_tufA.out 2>> screen_tufA.err

# Merge split libraries files
echo "Merging and converting back to fastq"
mkdir 06_Split_Library/tufA
cat 06_Split_Library/run1/seqs.fna 06_Split_Library/run2/seqs.fna 06_Split_Library/run3/seqs.fna > 06_Split_Library/tufA/seqs.fna
cat 06_Split_Library/run1/seqs_filtered.qual 06_Split_Library/run2/seqs_filtered.qual 06_Split_Library/run3/seqs_filtered.qual > 06_Split_Library/tufA/seqs_filtered.qual

# Merge qual and fasta again (Qiime)
convert_fastaqual_fastq.py -f 06_Split_Library/tufA/seqs.fna -q 06_Split_Library/tufA/seqs_filtered.qual -c fastaqual_to_fastq -o 06_Split_Library/tufA/

# Merge mapping file
merge_mapping_files.py -m 05_Mapping/R1_Map_tufA.txt,05_Mapping/R2_Map_tufA.txt,05_Mapping/R3_Map_tufA.txt -o 05_Mapping/Map_tufA.txt

# Quality filter using USEARCH (also trims at specific length with -fastq_trunclen, here 250 for tufA)
echo "Quality filter - Usearch"
mkdir 07_QC
usearch8 -fastq_filter 06_Split_Library/tufA/seqs.fastq -fastq_maxee 0.5 -fastqout 07_QC/1_all_recs_tufa_qc.fastq

# Dereplication
usearch8 -derep_fulllength 07_QC/1_all_recs_tufa_qc.fastq -fastqout 07_QC/2_no_replicates.fastq -sizeout

# Remove clusters with less than 5 reads
usearch8 -sortbysize 07_QC/2_no_replicates.fastq -fastqout 07_QC/3_all_good_seqs.fastq -minsize 5

# Pick OTUs with UPARSE (-otu_radius_pct 2 = ca. 98% similarity threshold)
echo "OTU clustering - Uparse"
mkdir 08_OTUs
usearch8 -cluster_otus 07_QC/3_all_good_seqs.fastq -relabel OTU_ -sizein -sizeout -otus 08_OTUs/otus.fasta -otu_radius_pct 2

usearch8 -usearch_global 07_QC/1_all_recs_tufa_qc.fastq -db 08_OTUs/otus.fasta -strand plus -id 0.98 -uc 08_OTUs/map.uc
	
# Assign taxonomy to the OTUs (QIIME)
echo "Assign taxomony - rdp - Qiime"
mkdir 09_Taxonomy
assign_taxonomy.py -i 08_OTUs/otus.fasta -t Ref_dataset_tufA/id_to_taxonomy_map.txt -r Ref_dataset_tufA/reference_sequences.fna -o 09_Taxonomy/rdp_assigned_taxonomy -m rdp
	
#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
echo "OTU table"
python bmp-map2qiime.py 08_OTUs/map.uc > 08_OTUs/otu_table.txt
make_otu_table.py -i 08_OTUs/otu_table.txt -t 09_Taxonomy/rdp_assigned_taxonomy/otus_tax_assignments.txt -o 08_OTUs/OTU_table.biom
biom summarize_table -i 08_OTUs/OTU_table.biom -o 08_OTUs/OTUs_stats1.txt


### Filter OTU table by min. count of reads per OTU per sample. = 5 
# Requires filter_observations_by_sample.py (https://gist.github.com/adamrp/7591573). Adam did a script to work with Qiime 1.9 as well.
echo "Filtering by sample"
filter_observations_by_sample.py -i 08_OTUs/OTU_table.biom -n 5 -o 08_OTUs/OTU_table_filter1.biom >> screen_tufA.out 2>> screen_tufA.err

# Filter by taxonomy: only o_Bryopsidales
filter_taxa_from_otu_table.py -i 08_OTUs/OTU_table_filter1.biom -p o_Bryopsidales -o 08_OTUs/OTU_table_filter2.biom >> screen_tufA.out 2>> screen_tufA.err

# Filter samples with less than 5 counts and controls:
filter_samples_from_otu_table.py -i 08_OTUs/OTU_table_filter2.biom -o 08_OTUs/OTU_table_filter3.biom -m 05_Mapping/Map_tufA.txt -s 'Host:*,!control' -n 5 >> screen_tufA.out 2>> screen_tufA.err


### Check with OTU stats and heatmap:
echo "Making stats and heatmaps"
biom summarize_table -i 08_OTUs/OTU_table_filter3.biom -o 08_OTUs/OTUs_stats_final.txt >> screen_tufA.out 2>> screen_tufA.err
#make_otu_heatmap_html.py -i 08_OTUs/OTU_table_filter3.biom -o 08_OTUs/OTU_Heatmap >> screen_tufA.out 2>> screen_tufA.err
make_otu_heatmap.py -i 08_OTUs/OTU_table_filter3.biom -o 08_OTUs/OTU_Heatmap >> screen_tufA.out 2>> screen_tufA.err


### Filter reference sequences:
echo "Filtering Reference sequences..."
filter_fasta.py -f 08_OTUs/otus.fasta -b 08_OTUs/OTU_table_filter3.biom -o 08_OTUs/rep_set_final.fna >> screen_tufA.out 2>> screen_tufA.err


echo "done"


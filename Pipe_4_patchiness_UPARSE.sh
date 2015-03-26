#!/bin/bash

# Pipeline for analysing environmental samples.
# Part 1 - preparing data for Qiime
# 26 - March - 2015

# Analyses of the patchiness experiment
# Run 2

echo "Launching....."


### Separate files by amplicon
echo "Separating files by amplicon and doing global trimming..."
### Separate genes (need separate_genes_w_global_trimming.py)
for file in 01_Extended_Frags/*.fastq; do separate_genes_w_global_trimming.py $file; done >> screen.out 2>> screen.err

mkdir 02_Separated_genes
mkdir 02_Separated_genes/16S
mkdir 02_Separated_genes/18S
mkdir 02_Separated_genes/UPA
mkdir 02_Separated_genes/tufA
mkdir 02_Separated_genes/Unmatched

mv 01_Extended_Frags/*.fastq_16S.fastq 02_Separated_genes/16S
mv 01_Extended_Frags/*.fastq_18S.fastq 02_Separated_genes/18S
mv 01_Extended_Frags/*.fastq_UPA.fastq 02_Separated_genes/UPA
mv 01_Extended_Frags/*.fastq_tufA.fastq 02_Separated_genes/tufA
mv 01_Extended_Frags/*.fastq_unmatched.fastq 02_Separated_genes/Unmatched

### Convert fastq to fasta and quality files (qiime) (or macqiime):
echo "Converting fastq to fasta and quality files..."
mkdir 03_Fasta_Qual
mkdir 03_Fasta_Qual/16S
mkdir 03_Fasta_Qual/18S
mkdir 03_Fasta_Qual/UPA
mkdir 03_Fasta_Qual/tufA

#only need if running macqiime
source /macqiime/configs/bash_profile.txt

for file in 02_Separated_genes/16S/*.fastq; do convert_fastaqual_fastq.py -f $file -c fastq_to_fastaqual -o 03_Fasta_Qual/16S/; done >>screen.out 2>> screen.err
for file in 02_Separated_genes/18S/*.fastq; do convert_fastaqual_fastq.py -f $file -c fastq_to_fastaqual -o 03_Fasta_Qual/18S/; done >>screen.out 2>> screen.err
for file in 02_Separated_genes/UPA/*.fastq; do convert_fastaqual_fastq.py -f $file -c fastq_to_fastaqual -o 03_Fasta_Qual/UPA/; done >>screen.out 2>> screen.err
for file in 02_Separated_genes/tufA/*.fastq; do convert_fastaqual_fastq.py -f $file -c fastq_to_fastaqual -o 03_Fasta_Qual/tufA/; done >>screen.out 2>> screen.err

# move unmatched files to a separate folder
mkdir 03_Fasta_Qual/16S/Undertemined_gene
mv 03_Fasta_Qual/16S/Undetermined* 03_Fasta_Qual/16S/Undertemined_gene
mkdir 03_Fasta_Qual/18S/Undertemined_gene
mv 03_Fasta_Qual/18S/Undetermined* 03_Fasta_Qual/18S/Undertemined_gene
mkdir 03_Fasta_Qual/UPA/Undertemined_gene
mv 03_Fasta_Qual/UPA/Undetermined* 03_Fasta_Qual/UPA/Undertemined_gene
mkdir 03_Fasta_Qual/tufA/Undertemined_gene
mv 03_Fasta_Qual/tufA/Undetermined* 03_Fasta_Qual/tufA/Undertemined_gene

#only need if running macqiime, will return an error in Linux, doesn't matter.
#source ~/.bash_profile

### Merge fasta and qual - need samples_list.csv with indexes
echo "merging..."
mkdir 04_All_records

#16S
mkdir 04_All_records/16S
merge_fasta.py 00_samples_list.csv 03_Fasta_Qual/16S/*.fna
merge_qual.py 03_Fasta_Qual/16S/*.qual
mv all_records.* 04_All_records/16S
mv map.txt 04_All_records/16S

#18S
mkdir 04_All_records/18S
merge_fasta.py 00_samples_list.csv 03_Fasta_Qual/18S/*.fna
merge_qual.py 03_Fasta_Qual/18S/*.qual
mv all_records.* 04_All_records/18S
mv map.txt 04_All_records/18S

#UPA
mkdir 04_All_records/UPA
merge_fasta.py 00_samples_list.csv 03_Fasta_Qual/UPA/*.fna
merge_qual.py 03_Fasta_Qual/UPA/*.qual
mv all_records.* 04_All_records/UPA
mv map.txt 04_All_records/UPA

#tufA
mkdir 04_All_records/tufA
merge_fasta.py 00_samples_list.csv 03_Fasta_Qual/tufA/*.fna
merge_qual.py 03_Fasta_Qual/tufA/*.qual
mv all_records.* 04_All_records/tufA
mv map.txt 04_All_records/tufA


### Quality control using Split Library. It is used here mainly to add labels to the sequences, the quality control will be done using UPARSE.
# Increased minimum lenght to (ensure the global alignment thing, 50bp less than actual size)
# parameters: min length (-l), quality score (-s = 1), max-homopolymer = 10. Save qual file!!
echo "Split_libraries running"
split_libraries.py -m 05_Mapping/16S/Map_16S.txt -f 04_All_records/16S/all_records.fna -q 04_All_records/16S/all_records.qual -p -s 1 -l 200 -b 16 -H 10 -o 06_Split_Library/16S --record_qual_scores >> screen_16S.out 2>> screen_16S.err
split_libraries.py -m 05_Mapping/18S/Map_18S.txt -f 04_All_records/18S/all_records.fna -q 04_All_records/18S/all_records.qual -p -s 1 -l 300 -b 16 -H 10 -o 06_Split_Library/18S --record_qual_scores >> screen_18S.out 2>> screen_18S.err
split_libraries.py -m 05_Mapping/UPA/Map_UPA.txt -f 04_All_records/UPA/all_records.fna -q 04_All_records/UPA/all_records.qual -p -s 1 -l 300 -b 16 -H 10 -o 06_Split_Library/UPA --record_qual_scores >> screen_UPA.out 2>> screen_UPA.err
split_libraries.py -m 05_Mapping/tufA/Map_tufA.txt -f 04_All_records/tufA/all_records.fna -q 04_All_records/tufA/all_records.qual -p -s 1 -l 400 -b 16 -H 10 -o 06_Split_Library/tufA --record_qual_scores >> screen_tufA.out 2>> screen_tufA.err

echo "Converting to fastq"
# Merge qual and fasta again (Qiime)
convert_fastaqual_fastq.py -f 06_Split_Library/16S/seqs.fna -q 06_Split_Library/16S/seqs_filtered.qual -c fastaqual_to_fastq -o 06_Split_Library/16S >> screen_16S.out 2>> screen_16S.err
convert_fastaqual_fastq.py -f 06_Split_Library/18S/seqs.fna -q 06_Split_Library/18S/seqs_filtered.qual -c fastaqual_to_fastq -o 06_Split_Library/18S >> screen_18S.out 2>> screen_18S.err
convert_fastaqual_fastq.py -f 06_Split_Library/UPA/seqs.fna -q 06_Split_Library/UPA/seqs_filtered.qual -c fastaqual_to_fastq -o 06_Split_Library/UPA >> screen_UPA.out 2>> screen_UPA.err
convert_fastaqual_fastq.py -f 06_Split_Library/tufA/seqs.fna -q 06_Split_Library/tufA/seqs_filtered.qual -c fastaqual_to_fastq -o 06_Split_Library/tufA >> screen_tufA.out 2>> screen_tufA.err

# Quality filter using USEARCH (also trims at specific length with -fastq_trunclen, here 250 for tufA)
echo "Quality filter - Usearch"
mkdir 07_QC
mkdir 07_QC/16S
mkdir 07_QC/18S
mkdir 07_QC/UPA
mkdir 07_QC/tufA

usearch8 -fastq_filter 06_Split_Library/16S/seqs.fastq -fastq_maxee 0.5 -fastqout 07_QC/16S/1_all_recs_qc.fastq >> screen_16S.out 2>> screen_16S.err
usearch8 -fastq_filter 06_Split_Library/18S/seqs.fastq -fastq_maxee 0.5 -fastqout 07_QC/18S/1_all_recs_qc.fastq >> screen_18S.out 2>> screen_18S.err
usearch8 -fastq_filter 06_Split_Library/UPA/seqs.fastq -fastq_maxee 0.5 -fastqout 07_QC/UPA/1_all_recs_qc.fastq >> screen_UPA.out 2>> screen_UPA.err
usearch8 -fastq_filter 06_Split_Library/tufA/seqs.fastq -fastq_maxee 0.5 -fastqout 07_QC/tufA/1_all_recs_qc.fastq >> screen_tufA.out 2>> screen_tufA.err

# Dereplication
usearch8 -derep_fulllength 07_QC/16S/1_all_recs_qc.fastq -fastqout 07_QC/16S/2_no_replicates.fastq -sizeout >> screen_16S.out 2>> screen_16S.err
usearch8 -derep_fulllength 07_QC/18S/1_all_recs_qc.fastq -fastqout 07_QC/18S/2_no_replicates.fastq -sizeout >> screen_18S.out 2>> screen_18S.err
usearch8 -derep_fulllength 07_QC/UPA/1_all_recs_qc.fastq -fastqout 07_QC/UPA/2_no_replicates.fastq -sizeout >> screen_UPA.out 2>> screen_UPA.err
usearch8 -derep_fulllength 07_QC/tufA/1_all_recs_qc.fastq -fastqout 07_QC/tufA/2_no_replicates.fastq -sizeout >> screen_tufA.out 2>> screen_tufA.err


# Remove clusters with less than 5 reads
usearch8 -sortbysize 07_QC/16S/2_no_replicates.fastq -fastqout 07_QC/16S/3_all_good_seqs.fastq -minsize 5 >> screen_16S.out 2>> screen_16S.err
usearch8 -sortbysize 07_QC/16S/2_no_replicates.fastq -fastqout 07_QC/16S/3_all_good_seqs.fastq -minsize 5 >> screen_18S.out 2>> screen_18S.err
usearch8 -sortbysize 07_QC/UPA/2_no_replicates.fastq -fastqout 07_QC/UPA/3_all_good_seqs.fastq -minsize 5 >> screen_UPA.out 2>> screen_UPA.err
usearch8 -sortbysize 07_QC/tufA/2_no_replicates.fastq -fastqout 07_QC/tufA/3_all_good_seqs.fastq -minsize 5 >> screen_tufA.out 2>> screen_tufA.err


# Pick OTUs with UPARSE (-otu_radius_pct 2 = ca. 98% similarity threshold for tufA. All other genes - default = 97%)
echo "OTU clustering - Uparse"
mkdir 08_OTUs
mkdir 08_OTUs/16S
mkdir 08_OTUs/18S
mkdir 08_OTUs/UPA
mkdir 08_OTUs/tufA

usearch8 -cluster_otus 07_QC/16S/3_all_good_seqs.fastq -relabel OTU_ -sizein -sizeout -otus 08_OTUs/16S/otus.fasta >> screen_16S.out 2>> screen_16S.err
usearch8 -cluster_otus 07_QC/18S/3_all_good_seqs.fastq -relabel OTU_ -sizein -sizeout -otus 08_OTUs/18S/otus.fasta >> screen_18S.out 2>> screen_18S.err
usearch8 -cluster_otus 07_QC/UPA/3_all_good_seqs.fastq -relabel OTU_ -sizein -sizeout -otus 08_OTUs/UPA/otus.fasta >> screen_UPA.out 2>> screen_UPA.err
usearch8 -cluster_otus 07_QC/tufA/3_all_good_seqs.fastq -relabel OTU_ -sizein -sizeout -otus 08_OTUs/tufA/otus.fasta -otu_radius_pct 2 >> screen_tufA.out 2>> screen_tufA.err

usearch8 -usearch_global 07_QC/16S/1_all_recs_qc.fastq -db 08_OTUs/16S/otus.fasta -strand plus -id 0.97 -uc 08_OTUs/16S/map.uc >> screen_16S.out 2>> screen_16S.err
usearch8 -usearch_global 07_QC/18S/1_all_recs_qc.fastq -db 08_OTUs/18S/otus.fasta -strand plus -id 0.97 -uc 08_OTUs/18S/map.uc >> screen_18S.out 2>> screen_18S.err
usearch8 -usearch_global 07_QC/UPA/1_all_recs_qc.fastq -db 08_OTUs/UPA/otus.fasta -strand plus -id 0.97 -uc 08_OTUs/UPA/map.uc >> screen_UPA.out 2>> screen_UPA.err
usearch8 -usearch_global 07_QC/tufA/1_all_recs_qc.fastq -db 08_OTUs/tufA/otus.fasta -strand plus -id 0.98 -uc 08_OTUs/tufA/map.uc >> screen_tufA.out 2>> screen_tufA.err
	
# Assign taxonomy to the OTUs (QIIME)
echo "Assign taxomony - rdp and uclust - Qiime"
mkdir 09_Taxonomy
mkdir 09_Taxonomy/16S
mkdir 09_Taxonomy/18S
mkdir 09_Taxonomy/UPA
mkdir 09_Taxonomy/tufA

assign_taxonomy.py -i 08_OTUs/16S/otus.fasta -t gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -r gg_13_8_otus/rep_set/97_otus.fasta -o 09_Taxonomy/16S/uclust_assigned_taxonomy >> screen_16S.out 2>> screen_16S.err
assign_taxonomy.py -i 08_OTUs/18S/otus.fasta -t Silva_111/taxonomy/97_Silva_111_taxa_map.txt -r Silva_111/rep_set/97_Silva_111_rep_set.fasta -o 09_Taxonomy/18S/uclust_assigned_taxonomy >> screen_18S.out 2>> screen_18S.err
assign_taxonomy.py -i 08_OTUs/UPA/otus.fasta -t Ref_dataset_UPA/id_to_taxonomy_map.txt -r Ref_dataset_UPA/reference_sequences.fna -o 09_Taxonomy/UPA/rdp_assigned_taxonomy -m rdp >> screen_UPA.out 2>> screen_UPA.err
assign_taxonomy.py -i 08_OTUs/tufA/otus.fasta -t Ref_dataset_tufA/id_to_taxonomy_map.txt -r Ref_dataset_tufA/reference_sequences.fna -o 09_Taxonomy/tufA/rdp_assigned_taxonomy -m rdp >> screen_tufA.out 2>> screen_tufA.err


#align_seqs.py -i 08_OTUs/18S/rep_set.fna -t Silva_111/rep_set_aligned/97_Silva_111_rep_set.fasta -o 09_Alignment/18S/ >> screen_18S.out 2>> screen_18S.err
#filter_alignment.py -i 09_Alignment/18S/rep_set_aligned.fasta -s -o 09_Alignment/18S/ >> screen_18S.out 2>> screen_18S.err

	
#Convert UC to otu-table.txt <<< BMP SCRIPT>>>
#echo "OTU table"
#python bmp-map2qiime.py 08_OTUs/map.uc > 08_OTUs/otu_table.txt
#make_otu_table.py -i 08_OTUs/otu_table.txt -t 09_Taxonomy/rdp_assigned_taxonomy/otus_tax_assignments.txt -o 08_OTUs/OTU_table.biom
#biom summarize_table -i 08_OTUs/OTU_table.biom -o 08_OTUs/OTUs_stats1.txt






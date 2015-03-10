#!/bin/bash

# Pipeline for analysing environmental samples.
# Part 1 - preparing data for Qiime
# 19 - Jan - 2015

echo "Launching....."


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

echo "End of pipe 1"
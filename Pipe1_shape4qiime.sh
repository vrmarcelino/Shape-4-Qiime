#!/bin/bash

# Pipeline for analysing environmental samples.
# Part 1 - preparing data for Qiime
# 19 - Jan - 2015

echo "Launching....."

### Trim adaptors (need trim_adaptors.py)
#Check the forward and reverse file names and change accordingly if needed.
for file in *_R1_001.fastq; do trim_adaptors.py $file forward; done >> screen_for.out 2>> screen_for.err
for file in *_R2_001.fastq; do trim_adaptors.py $file reverse; done >> screen_rev.out 2>> screen_rev.err


### Merge paired-end files using flash
# Options: -m = minimum overlap / -M = maximum overlap / -d = output directory / -O because of the adaptor trimmings / -x = ratio of allowed mismatches
# $flash VRM135-M3a_S8_L001_R1_001.fastq_trimmed.fastq VRM135-M3a_S8_L001_R2_001.fastq_trimmed.fastq -m 70 -M 305 -O -x 0.1 -d Merged -o _paired
# for all of them, run Heroen's perl (need to edit file names if different!)
run_flash.pl

### Move extended frags (output from FLASH) to a folder and delete trash
mkdir 01_Extended_Frags
mv Merged/*extendedFrags.fastq 01_Extended_Frags
rm -r Merged
rm *fastq_trimmed.fastq

### Separate genes (need separate_genes.py)
for file in 01_Extended_Frags/*.fastq; do separate_genes.py $file; done >> screen.out 2>> screen.err

### Separate files by amplicon
echo "Separating files by amplicon..."
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

#only need if running macqiime, will return an error in Linux, doesn't matter.
source ~/.bash_profile

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
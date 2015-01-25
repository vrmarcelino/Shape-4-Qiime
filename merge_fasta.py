# -*- coding: utf-8 -*-
""" Concatenate different fasta files and add barcodes.
Run this script after separate fasta and qual files (see onvert_fastaqual_fastq.py from qiime)

Created on Thu Jul 31 15:49:39 2014

@author: VanessaRM

Still need to be done: match the .fna sample name with the sample_ID in the csv file.
At this point, this script will add the indexes by alphabetic(?) order, 
so the indexe-sample match is not the same as the orginal ones.


"""

from Bio import SeqIO
import sys
import pandas # for .csv handling

#input files:
samples_indexes = str(sys.argv[1])
si = pandas.read_csv("samples_list_run2.csv")

Sample_ID = si["Sample_ID"]
Index_ID = si["Frd_index_name"] + si["Rev_index_name"]
Index_seq = si["Frd_Index"] + si["Rev_Index_RC"]


input_fasta = []
for n in sys.argv[1:]:
    input_fasta.append(str(n))

#Store the files
all_records = []

#function for adding barcode sequences
def add_barcode(records, barcode):
    for seq_record in records:
        seq_record.seq = (barcode + seq_record.seq)
        all_records.append (seq_record)

#iterate over input files
counter = 0
for file in input_fasta:
    original_reads = SeqIO.parse(file, "fasta")

# match sample name
    full_sample_name = str(file)
    name_split = full_sample_name.split("_")
    sample_name = name_split[0]
    
    barcode_seq = Index_seq[counter]
    print""
    print "Adding the barcode %s to the %s file" %(barcode_seq, file)
    do_it = add_barcode(original_reads, barcode_seq)
    counter +=1
  

SeqIO.write(all_records, "all_records.fna", "fasta")
print""
print "Done!"

    
        
        
        
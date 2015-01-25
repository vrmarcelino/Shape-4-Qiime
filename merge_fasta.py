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
si = pandas.read_csv(samples_indexes)
Index_seq = si["Frd_Index"] + si["Rev_Index_RC"]

input_fasta = []
for n in sys.argv[2:]:
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
mapping_file = ["#SampleID"+'\t'+"BarcodeSequence"+'\t'+"LinkerPrimerSequence"+'\t'+"BlaBlaBla"+'\t'+"Description"]
for file in input_fasta:
    original_reads = SeqIO.parse(file, "fasta")

    barcode_seq = Index_seq[counter]
    print""
    print "Adding the barcode %s to the %s file" %(barcode_seq, file)
    do_it = add_barcode(original_reads, barcode_seq)
    
# Store info for mapping file
    full_sample_name = str(file)
    name_split = full_sample_name.split("_")
    sample_name = name_split[0]
    mapping_file.append(sample_name + '\t' + barcode_seq)
     
    counter +=1

# Save stuff
SeqIO.write(all_records, "all_records.fna", "fasta")

savefile = open("map.txt", "w")
for lines in mapping_file:
    savefile.write("%s\n" % lines)

print""
print "Mapping file saved as 'map.txt'"

print ""

print "Done!"

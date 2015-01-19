# -*- coding: utf-8 -*-
"""
Demultiplexing samples based on index sequences
Allows 1 mismatch for each index.

Created on Tue Aug  5 10:58:40 2014
@author: VanessaRM
"""
from Bio import SeqIO
import regex
import sys
import pandas
from itertools import izip

#Help text
if len(sys.argv) < 5:
    print ""
    print "Script to demultiplex paired end fastQ files."
    print "" 
    print "Usage: supply the forward and reverse files, indexes, a .csv file with samples and indexes, and the maximum allowed index mismatch"
    print ""    
    print "ex: python demultiplex.py Frw_reads_L001_R1_001.fastq Rrvr_eads_L001_R2_001.fastq Frd_indexes_ L001_I1_001.fastq Rvr_indexes_ L001_I2_001.fastq samples.csv 1"
    print ""
    sys.exit()


#input files (note that the L001 matches with I2 because the Illumina adaptors have been inversed!)
input_file_R1 = str(sys.argv[1])
input_file_R2 = str(sys.argv[2])

input_index_I1 = str(sys.argv[3]) 
input_index_I2 = str(sys.argv[4])

samples_indexes = str(sys.argv[5])

max_mismatch = str(sys.argv[6])



def classificator (seq_f, seq_r, samples_df, mismatches = 1):
    "This function takes one cluster (i.e. 2 index sequences - forward and reverse)"
    " and returns the sample bin to wich that index pair belongs to."
    
    for i in xrange(len(samples_df)):
    
        f_index = "(" + samples_f_ind[i] + ")" + "{e<=" + str(mismatches) + "}"
        search_f = regex.search((f_index), seq_f)
        
        r_index = "(" + samples_r_ind[i] + ")" + "{e<=" + str(mismatches) + "}"
        search_r = regex.search((r_index), seq_r)
        
        if search_f and search_r != None:
            return samples_id[i]


def demultiplex (f_reads, f_indexes, r_reads, r_indexes, samples_df, mismatches):
    "Function to demultiplex samples"
    for Fr, Fi, Rr, Ri in izip(f_reads, f_indexes, r_reads, r_indexes):
        
        # Match the indexes to a sample
        Fi_seq = str(Fi.seq) #forwrad index seq
        Ri_seq = str(Ri.seq) #Rev index seq
        classify = classificator(Fi_seq, Ri_seq, samples_df, mismatches)
        
        # Save reads
        dict_forward_reads[classify].append(Fr)
        dict_reverse_reads[classify].append(Rr)


# Read input files:
f_reads = SeqIO.parse(input_file_R1, "fastq")
f_indexes = SeqIO.parse(input_index_I1, "fastq")

r_reads = SeqIO.parse(input_file_R2, "fastq")
r_indexes = SeqIO.parse(input_index_I2, "fastq")

si = pandas.read_csv(samples_indexes)

samples_id = si["Samples"]
samples_f_ind = si["Frd_Index"]
samples_r_ind = si["Rev_Index_RC"]


# Create dictionaries containing lists to store the separated samples
dict_forward_reads = {}
for i in samples_id:
    dict_forward_reads['%s' %i] = []
    dict_forward_reads[None] = []
    
dict_reverse_reads = {}
for i in samples_id:
    dict_reverse_reads['%s' %i] = []
    dict_reverse_reads[None] = []

print ""
print "demultiplexing..."
do_it = demultiplex(f_reads, f_indexes, r_reads, r_indexes, si, max_mismatch)


#Save files
count_all = [] # check if all seqs are being saved
print ""
for sample in dict_forward_reads:   
    file_name_f = '%s' %sample + "_Forward.fastq"    
    count =  SeqIO.write(dict_forward_reads[sample], file_name_f, "fastq")
    count_all.append(count)
    print "Saved %i reads in the %s file." % (count, file_name_f)

print ""
for sample in dict_reverse_reads:   
    file_name_f = '%s' %sample + "_Reverse.fastq"    
    count =  SeqIO.write(dict_forward_reads[sample], file_name_f, "fastq")
    print "Saved %i reads in the %s file." % (count, file_name_f)

print ""
print "A total of %i reads have been saved." % sum(count_all)
print "Done!"

#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Trim tail and adaptors from seq sequences.
The 16S is a short amplicon, so the adaptor and tail has been sequenced in the final end of these sequences. 
We need therefore to trim them in order to merge the paired end reads.

Created on Jul 30 2014
modified: Jan 2015
"""

from Bio import SeqIO
import regex
import sys

# Help
if len(sys.argv) == 1:
    print ""
    print "Script to trimm the tail and adaptors from short reads"
    print ""
    print "Usage: supply the file name and inform if its forward or reverse"
    print "ex: python trim_adaptors.py my_fastq_file.fastq forward"
    print ""
    print "If you want to run it for multiple files, use the shell:"
    print "for file in *_R1_001.fastq; do python trim_adaptors.py $file forward; done >> screen.out 2>> screen.err &"
    print ""
    sys.exit()
    
    
# input files and arguments:
input_file = str(sys.argv[1])
output_file = input_file + "_trimmed.fastq"


#Tail in forward reads:
#full_tail_f = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" # Reverse complement of the reverse tail
tail_small_f = '(CTGTCTCTTATACAC){e<=1}' # first 15 bases, allowing 1 mismatch

#Tail in reverse reads:
#full_tail_r = "CTGTCTCTTATACACATCTGACGCTGCCGACGA" # reverse complement of the forward tail
tail_small_r = "(CTGTCTCTTATACAC){e<=1}" # first 15 bases, allowing 1 mismatch


#define trimmer function
def trimmer (records, tail):
    "Trims the tail and all the bases after it (adaptors and etc...)"
    for record in records:
        sequence = str(record.seq)        
#        index = record.seq.find(tail)
        index = regex.search((tail), sequence)       
        if index == None:
            #tail not found
            yield record
        else:
            # Trim it!
            cut_off = (int(index.span()[0]) - 6) #where the tail starts minus 6 bp to cut off linkers and spacers
            yield record [:cut_off]
        
        

original_reads = SeqIO.parse(input_file, "fastq")

if str(sys.argv[2]) == "forward":
    trimmed_reads = trimmer(original_reads, tail_small_f)
    
elif str(sys.argv[2]) == "reverse":
    trimmed_reads = trimmer(original_reads, tail_small_r)

#Else return error:
else:
    print ("")
    print ("Please define if the file is forward or reverse!!")
   

count = SeqIO.write(trimmed_reads, output_file, "fastq")
print""
print "Saved %i reads in the %s file." % (count, output_file)



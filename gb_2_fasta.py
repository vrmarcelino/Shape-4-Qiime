#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Converts genbank files to fasta
VRM
August 2015
'''

from Bio import SeqIO
import sys
 
input_reference_dataset = str(sys.argv[1])
output = str(str(input_reference_dataset) + ".fasta")

input_handle = open(input_reference_dataset, "rU")
output_handle = open(output, "w")
 
sequences = SeqIO.parse(input_handle, "genbank")
count = SeqIO.write(sequences, output_handle, "fasta")


output_handle.close()
input_handle.close()
print "Done! Converted %i records" % count



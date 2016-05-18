#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to find which sequences were filtered by trimAl
Procudes a fasta file with the filtered seqs.
Created on Fri Aug 21 10:43:58 2015

@author: VanessaRM
"""
from argparse import ArgumentParser
from Bio import SeqIO
import sys
#from Bio.SeqRecord import SeqRecord

#import os
#os.chdir('/Users/VanessaRM/Documents/PhD/4_Merged_123Runs/05_alignments')

parser = ArgumentParser()
parser.add_argument('-oa', '--original_aln', help='The path to the input aligned pfiltered OTUs', required=True)
parser.add_argument('-fa', '--filtered_aln', help='The path to the filtered OTUs - output from trimAL.txt', required=True)
parser.add_argument('-o', '--output', help='The path and name of the output file', required=True)

args = parser.parse_args()
original_aln = args.original_aln
filtered_aln = args.filtered_aln
output_fp = args.output

all_seqs_in_filtered_aln = []
store_bad_otus = []


#original_aln = "UPA.nt_aligned.fas.reduced_pfiltered.fasta"
#filtered_aln = "test_trimal_UPA_80_80.fasta"

# get a list of OTUs in filtered alignment:
for seq_record in SeqIO.parse(filtered_aln, "fasta"):
    all_seqs_in_filtered_aln.append(seq_record.id)

# Check which OTUs are not in the filtered aln
for seq_record in SeqIO.parse(original_aln, "fasta"):
    if seq_record.id not in all_seqs_in_filtered_aln:
        store_bad_otus.append(seq_record)


SeqIO.write(store_bad_otus, str(output_fp), "fasta")
print ""
print "%i sequences saved in %s" %(len(output_fp), output_fp)


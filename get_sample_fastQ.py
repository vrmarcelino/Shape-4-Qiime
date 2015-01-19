# -*- coding: utf-8 -*-
"""
Get a few sequences from a huge fastQ file to make some tests
Created on Tue Jul 15 10:56:12 2014
@author: VanessaRM
"""

from Bio import SeqIO

wish_list = 500
sample_list = []
min_seq_lengh = 100 #bp

for rec in SeqIO.parse("Undetermined_S0_L001_R1_001.fastq", "fastq"):    
    if len(rec.seq) > min_seq_lengh:
        if len(sample_list) < wish_list:
            sample_list.append(rec)
        else:
            break


reads = SeqIO.write(sample_list, "sample_fastQ_500sq.fastq", "fastq")
print "Saved %i reads" % reads


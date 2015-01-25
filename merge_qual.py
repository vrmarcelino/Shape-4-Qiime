#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Concatenate the quality scores files and add a quality for the barcodes

Created on Thu Jul 31 18:38:36 2014

@author: VanessaRM
"""

from Bio import SeqIO
import re
import sys


all_qual =[]


#input files:
input_files = []
for n in sys.argv[1:]:
    input_files.append(str(n))
    
    
barcode_qual = "30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 " # 16 good quality scores

def add_quality(records, barcode_qual):
    "Function add the quality scores to the barcode"
    for seq_record in records:
        original_qual = str(seq_record.letter_annotations["phred_quality"])
        new_qual = (barcode_qual + original_qual)
        new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
        name = ">" + (seq_record.description)
        
        all_qual.append (name)
        all_qual.append (new_qual)


#Iterate over input files:
for file in input_files:
    print""
    print "Processing %s file" %file
    original_reads = SeqIO.parse(file, "qual")
    do_it = add_quality(original_reads, barcode_qual)


#Save:
savefile = open("all_records.qual", "w")
for lines in all_qual:
    savefile.write("%s\n" % lines)


print "Done!"


#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Note - this still doesn't work properly. It writes the fasta file but stop looping throough the genbank file
Converts genbank files to fasta skipping sequences badly formatted in genbank
VRM
August 2015
'''


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import sys

#help
if len(sys.argv) == 1:
    print ""
    print "Takes a genbank dataset file and produces a fasta, ignoring entries that try to screw up everything"
    print ""
    print "Usage: gb_2_fasta_ignoring_bad_entries. py 23_full_dataset.gb"
    print ""
    print ""
    sys.exit()


input_reference_dataset = str(sys.argv[1])
output = str(str(input_reference_dataset) + ".fasta")
store_seq_rec = []

### Save record with ID and sequence only
try:
    for seq_record in SeqIO.parse(input_reference_dataset, "genbank"):
        seq = str(seq_record.seq)
        newrec=SeqRecord(Seq(seq),id=(seq_record.name),name="",description="")
        store_seq_rec.append(newrec)
except IndexError:
    print 'skipping bad record'

# Save file
count = SeqIO.write(store_seq_rec, output, "fasta")
print "Done! Converted %i records" % count


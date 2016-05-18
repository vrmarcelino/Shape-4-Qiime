# -*- coding: utf-8 -*-
"""
Filter taxonomic assignments file based on fasta file.
@author: VanessaRM
"""

from Bio import SeqIO
from argparse import ArgumentParser


parser = ArgumentParser()
parser.add_argument('-i', '--otus_fasta', help='The path to the input fasta file with OTU sequences', required=True)
parser.add_argument('-t', '--taxonomic_assignments', help='The path to the original taxonomic assignmnets', required=True)
parser.add_argument('-o', '--output', help='The path and name of the output file', required=True)

args = parser.parse_args()
otus_fasta = args.otus_fasta
tax_file = args.taxonomic_assignments
output = args.output
store_seqs_ids = []
wanted_taxons = []

for seq_record in SeqIO.parse(otus_fasta, "fasta"):
    full_id = seq_record.id
    store_seqs_ids.append(full_id)

tax = open(tax_file, "r" )
for line in tax:
    identifier = line.split('\t')
    otu = identifier[0]
    if otu in store_seqs_ids:
        wanted_taxons.append (line)


savefile = open(output, "w")
for i in wanted_taxons:
    savefile.write("%s\n" % i)

print ""
print "Done!"
print ""
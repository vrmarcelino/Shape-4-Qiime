# -*- coding: utf-8 -*-
"""
Exclude gene bank numbers from sequence identifiers.
Produces txt file with IDs
Created on Fri Aug 28 17:28:40 2015
@author: VanessaRM
"""


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

alignment = str(sys.argv[1])
output_aln = str(str(alignment) + ".new_ids.fas")
output_ids = str(str(alignment) + ".genbank_ids.txt")

new_aln = []
store_gb_ids = []

#filter gb ids and "size info from OTUs"
for seq_record in SeqIO.parse(alignment, "fasta"):
    full_id = seq_record.id
    store_gb_ids.append(seq_record.id)
    
    pieces = full_id.split('_')
    new_id = pieces[0] + "_" + pieces[1]
    seq = str(seq_record.seq)
    
    new_record = SeqRecord(Seq(seq), id= new_id, description='')
    
    # OTUs
    new_rec_2 = new_record.id
    pieces2 = new_rec_2.split(';')
    new_id_2 = pieces2[0]
    
    new_record_2 = SeqRecord(Seq(seq), id= new_id_2, description='')
    
    
    new_aln.append(new_record_2)
    
    
    
savefile = open(output_ids, "w")
for lines in store_gb_ids:
    savefile.write("%s\n" % lines)


count = SeqIO.write(new_aln, output_aln, "fasta")
print""
print "Done. Saved %i sequences in the %s file." % (count, output_aln)
print""


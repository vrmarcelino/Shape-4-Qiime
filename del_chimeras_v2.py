'''
Delete chimeras from OTU txt file for build OTU table in QIIME
Version 2 - use only the filtered seqs that you want to keep.
Created on 18/12/2014
input = Otus_filtered.fasta (from geneious)
		seqs_otus.txt (from 03_OTUs)

'''
from Bio import SeqIO

not_trash_id =[]
non_chimeric_OTUS = []

for seq_record in SeqIO.parse("Otus_filtered.fasta" , "fasta"):
    not_trash_id.append (seq_record.id)

ins = open( "seqs_otus.txt", "r" )
for line in ins:
    identifier = line.split('\t')
    denovo = identifier[0]
    if denovo in not_trash_id:
        non_chimeric_OTUS.append (line)
        
print (len(non_chimeric_OTUS))


savefile = open("clean_OTUS.txt", "w")
for lines in non_chimeric_OTUS:
    savefile.write("%s" % lines)


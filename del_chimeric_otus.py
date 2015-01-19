'''
Script to delete unwanted OTUS from both OTU_rep.fna and the seqs_otus.txt
Used here to get only Chlorophytes reads
Run this script after assign taxonomy with RDP.
Created on 19/08/2013
Modified on Aug 29
'''
#import
from Bio import SeqIO

#create lists
ulvo_denovo_list =[]
others_denovo_list =[]

ulvo_seqs = []
ulvo_otus =[]

#read RDP output
RDP_out = open("OTU_rep_tax_assignments.txt", "r")

#get the ulvophyceae IDs:
for line in RDP_out:
    get_denovo = line.split('\t')
    denovo_x = get_denovo[0]
    classification = get_denovo[1]

    if classification.find("Chlorophyta") != -1:
        ulvo_denovo_list.append(denovo_x)


#read the sequence file
for seq_record in SeqIO.parse("OTU_rep.fna" , "fasta"):
    if (seq_record.id) in ulvo_denovo_list:
        ulvo_seqs.append(seq_record)


print ("found %i Chlorophyta sequences" %len(ulvo_seqs))


#Delete otus from otu mapping file
otus_all = open( "seqs_otus.txt", "r" )
for line in otus_all:
    identifier = line.split('\t')
    denovo_y = identifier[0]
    if denovo_y in ulvo_denovo_list:
        ulvo_otus.append (line)

#save files
SeqIO.write(ulvo_seqs, "Chlorophytes_seqs.fna", "fasta")

savefile = open("Chlorophytes_otus.txt", "w")
for lines in ulvo_otus:
    savefile.write("%s" % lines)

print ("Done! Files saved as 'Chlorophytes_seqs.fna' and 'Chlorophytes_otus.txt'.")
'''
Script to delete unwanted OTUS from both OTU_rep.fna and the seqs_otus.txt
Used here to get only Ulvophyceae otus for the Porites rus analysis.
Run this script after align (and exlude) otus in Geneious and assign taxonomy with RDP.
obs: check for _(reversed) in the seqs names!
Created on 19/08/2013
'''
#import
from Bio import SeqIO

#create lists
ulvo_denovo_list =[]
others_denovo_list =[]

ulvo_seqs = []
ulvo_otus =[]

#read RDP output
RDP_out = open("Porites_2746_otus_tax_assignments.txt", "r")

#get the ulvophyceae IDs:
for line in RDP_out:
    get_denovo = line.split('\t')
    denovo_x = get_denovo[0]
    classification = get_denovo[1]

    if classification.find("Ulvophyceae") != -1:
        ulvo_denovo_list.append(denovo_x)


#read the sequence file
for seq_record in SeqIO.parse("OTU_rep.fna" , "fasta"):
    if (seq_record.id) in ulvo_denovo_list:
        ulvo_seqs.append(seq_record)


print ("found %i ulvophyceae sequences" %len(ulvo_seqs))


#Delete otus from otu mapping file
otus_all = open( "Porites_seqs_otus.txt", "r" )
for line in otus_all:
    identifier = line.split('\t')
    denovo_y = identifier[0]
    if denovo_y in ulvo_denovo_list:
        ulvo_otus.append (line)

#save files
SeqIO.write(ulvo_seqs, "Porites_ulvo_seqs.fna", "fasta")

savefile = open("Porites_ulvo_otus.txt", "w")
for lines in ulvo_otus:
    savefile.write("%s" % lines)

print ("Done! Files saved as 'Porites_ulvo_seqs.fna' and 'Porites_ulvo_otus.txt'.")
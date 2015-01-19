'''
Script to refine the reference datasets used to identify OTUs
Created on 23/07/2013
'''

# Check if there are repetitive seqs:
from Bio import SeqIO

all_species = []
unique_species =[]

for seq_record in SeqIO.parse("tufA_ulvophyceae_reference.gb" , "genbank"):
    all_species.append(seq_record.annotations["organism"])
    
for seq_record in SeqIO.parse("tufA_ulvophyceae_reference.gb" , "genbank"):
    if str(seq_record.annotations["organism"]) not in unique_species:
        unique_species.append(seq_record.annotations["organism"])

print ("number of all species in dataset is %i" %(len(all_species)))
print ("number of unique species is %i" %(len(unique_species)))



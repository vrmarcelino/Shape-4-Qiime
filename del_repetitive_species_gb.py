'''
Script to refine the reference datasets
Deletes the repetitive species in a dataset and save the unique sequence records in a .gb file
Created on 23/07/2013
'''

from Bio import SeqIO

all_species = []
unique_species = []
unique_records = []

# All species
for seq_record in SeqIO.parse("tufA_ulvophyceae_reference.gb" , "genbank"):
    all_species.append(seq_record.annotations["organism"])

# Find the unique species
for seq_record in SeqIO.parse("tufA_ulvophyceae_reference.gb" , "genbank"):
    if str(seq_record.annotations["organism"]) not in unique_species:
        unique_species.append(seq_record.annotations["organism"])
        unique_records.append(seq_record)

#print results
deleted_records = (len(all_species) - len(unique_records))
print ("number of deleted records is %i" %(deleted_records))
print ("%i species remained in the dataset" %(len(unique_records)))

# write to file
SeqIO.write(unique_records, "unique_records_tufA.gb", "genbank")

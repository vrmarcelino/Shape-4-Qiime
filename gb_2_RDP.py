'''
Files needed for RDP classifier:
Create files 2 files from a genbank reference dataset:
RDP requires a fASTA file of reference sequences and a taxonomy file that maps IDs of the reference seqs to a taxonomic hierarchy
eg. file 1 (reference_sequences_RDP): 
KC191578 aaagaacata ttttact....
file 2 (id_to_taxonomy_map)
HE600176    d_Eukaryota;k_Viridiplantae;p_Chlorophyta;c_Ulvophyceae;o_Ulvales;f_Ulvaceae;g_Ulva;s_Ulva_californica

WARNING: The classification (into domain/kingdom/class/etc) depends of the order and number of taxonomic ranks defined in the genbank record. 

Created on 15/08/2013
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
#import sys

# 1 ### Save file 1 (fasta with ID and sequence)
file1 = []

#for user input * sys.argc[1]
for seq_record in SeqIO.parse("UPA_reference_dataset.gb" , "genbank"):
    seq = str(seq_record.seq)
    newrec=SeqRecord(Seq(seq),id=(seq_record.name),name="",description="")
    file1.append(newrec)

#save
SeqIO.write(file1, "reference_sequences_RDP.fna", "fasta")


# 2 ### Make file 2 (id_to_taxonomy_map)
file2 = []

for seq_record in SeqIO.parse("UPA_reference_dataset.gb" , "genbank"):
    iden = str(seq_record.name)
    
    descr_raw = str(seq_record.annotations["taxonomy"])
    parts_raw = descr_raw.split(",")
    genera = parts_raw[-1]
    if len(parts_raw) < 6:
        descr_raw = str(seq_record.annotations["taxonomy"]) + (",x" * (6-(descr_raw.count(",")))) # Ensure that all records have at least 6 tax. levels 
    parts = descr_raw.split(",")
    
# Rhodophyta
    if parts[1].find ("'Rhodophyta'") != -1:
        descr = ("d_" + parts[0]+";")+("k_Plantae;")+("p_" + parts[1]+";")+("c_"+ parts[2]+";")+("o_"+parts[3]+";")+("f_"+parts[4]+";")+("g_"+ genera+";")

# Bacteria
    elif parts[0].find ("'Bacteria'") != -1:
        descr = ("d_Procaryota;")+("k_" + parts[0]+";")+("p_" + parts[1]+";")+("c_x;")+("o_"+parts[2]+";")+("f_"+parts[3]+";")+("g_"+ genera+";")

# Green algae and etc
    else:
        descr = ("d_" + parts[0]+";")+("k_" + parts[1]+";")+("p_" + parts[2]+";")+("c_"+ parts[3]+";")+("o_"+parts[4]+";")+("f_"+parts[5]+";")+("g_"+ genera+";")

    sp = str(seq_record.annotations["organism"])
    sp = re.sub(r" ","_",sp)
 
    file2.append(iden + '\t' + descr + "s_" + sp)


# Format the file - delete "[", "]" and replace "," for ";"
file2_2= []
    
for lines in file2:
    file2_1 = re.sub(r"\['|\.|\-|\'|\ |\]","",lines) # delete [' and dots and spaces and etc
    file2_2.append(file2_1)

#save
savefile = open("id_to_taxonomy_map.txt", "w")
for lines in file2_2:
    savefile.write("%s\n" % lines)

print ('Done! Files were saved as "reference_sequences_RDP.fna" and "id_to_taxonomy_map.txt".')

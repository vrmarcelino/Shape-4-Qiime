'''
Delete chimeras from OTU txt file for build OTU table in QIIME
The chimera_seqs.fna can be saved form the geneious map_to_reference (select "save report", click on the seqs that did not align and then export as fasta)
-remember to take out also the ones that aligned but were excluded (badly aligned)
Created on 13/08/2013

'''
from Bio import SeqIO

chimeras_id =[]
non_chimeric_OTUS = []

for seq_record in SeqIO.parse("chimeras_seqs.fas" , "fasta"):
    chimeras_id.append (seq_record.id)

ins = open( "tufA_seqs_otus.txt", "r" )
for line in ins:
    identifier = line.split('\t')
    denovo = identifier[0]
    if denovo not in chimeras_id:
        non_chimeric_OTUS.append (line)
        
print (len(non_chimeric_OTUS))


savefile = open("non_chimeric_OTUS.txt", "w")
for lines in non_chimeric_OTUS:
    savefile.write("%s" % lines)


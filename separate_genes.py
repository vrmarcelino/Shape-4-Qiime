separate_genes.py
# -*- coding: utf-8 -*-
"""
Script to separate genes from MiSeq run fastQ files (NEXTERA style)
AFTER the reads have been trimmed for the tail and adapter (see trim_adaptors.py)
Trim forward and reverse primers.
It will produce one file for each gene (4 here - 16S, 18S, UPA, tufA)

Created on Mon Jul 28 17:27:33 2014
@author: VanessaRM
"""
""" USAGE:
python speratae_genes.py my_fastq_file.fastq 
"""


from Bio import SeqIO
import regex
import sys

# Help text
if len(sys.argv) == 1:
    print ""
    print "Script to separate genes and trim primers from a fastQ file."
    print "Allows one mismatch in the primer sequence"
    print ""
    print "Usage: supply the file name"
    print "ex: python speratae_genes.py my_fastq_file.fastq"
    print ""
    print "If you want to run it for multiple files, use the shell:"
    print "for file in *_R1_001.fastq; do python separate_genes.py $file; done >> screen.out 2>> screen.err &"
    print ""
    sys.exit()



# input files and arguments:
input_file =  str(sys.argv[1])
output_file_16S = input_file + "_16S.fastq"
output_file_18S = input_file + "_18S.fastq"
output_file_UPA = input_file + "_UPA.fastq"
output_file_tufA = input_file + "_tufA.fastq"
output_file_umatched = input_file + "_unmatched.fastq"

# Globals
finder = ""

# Create lists to store seqs
#all_seqs =[] #just to see if the code behaves well

s16_list =[]
s18_list = []
UPA_list=[]
tufA_list =[]
unmatched=[]

# Forward Primers (\w = ambiguous positions, only 12 "internal" bases used), allowing for 1 mismatch.
p16S_f = '(GTGCCAGC\wGCCGCGGTAA){e<=1}' #F515
p18S_f = '(GGTGGTGCATGGCCGTTCTTAGTT){e<=1}' #NF1
pUPA_f = '(GGACAGAAAGACCCTATGAA){e<=1}' #p23SrV_f1
ptufA_f = '(AC\wGG\wCG\wGG\wAC\wGT){e<=1}' #Oq_tuf2F

#Reverse primers(their reverse complement):
p16S_r = '(ATTAGA\wACCC\w\wGTAGTCC){e<=1}' #R806 rc ATTAGAWACCCBDGTAGTCC
p18S_r = '(ATTACGTCCCTGCCCTTTGTA){e<=1}' # B18Sr2B rc
pUPA_r = '(CTCTAGGGATAACAGGCTGA){e<=1}' #p23SrV_r1 rc
ptufA_r = '(GCG\wTT\wGC\wATTCG\wGAAGG){e<=1}' #tufAR rc


#define primer_finder function
def primer_finder (records, primer_f1, primer_r1, primer_f2, primer_r2, primer_f3, primer_r3, primer_f4, primer_r4):
    "Trims the primers and saves the sequences in amplicon-separated files. Put primers in order: 16S, 18S, UPA, tufA"
    for record in records:
        sequence = str(record.seq)
        
        #Initial values of cut_off - in case you don't need to trim anything
        cut_off_f= 0
        cut_off_r = len(sequence)
        
        #Search the primers f and r
        index_f = regex.search((primer_f1), sequence)
        index_r = regex.search((primer_r1), sequence)
        if index_f != None:
            #found the forward primer, so define where the sequence needs to be trimmed
            cut_off_f = int(index_f.span()[1])
        if index_r != None:
            #found the reverse primer
            cut_off_r = int(index_r.span()[0]+1) #the +1 is to cut just when the primer starts.
        
        #Store the trimmed seq
        if index_f or index_r != None:
            s16_list.append(record [cut_off_f:cut_off_r])
            
        else: #search for next primer
            index_f = regex.search((primer_f2), sequence)
            index_r = regex.search((primer_r2), sequence)
            if index_f != None:
                cut_off_f = int(index_f.span()[1])
            if index_r != None:
                cut_off_r = int(index_r.span()[0]+1)
            if index_f or index_r != None:
                s18_list.append(record [cut_off_f:cut_off_r])
           
            else:
                index_f = regex.search((primer_f3), sequence)
                index_r = regex.search((primer_r3), sequence)
                if index_f != None:
                    cut_off_f = int(index_f.span()[1])
                if index_r != None:
                    cut_off_r = int(index_r.span()[0]+1)
                if index_f or index_r != None:
                    UPA_list.append(record [cut_off_f:cut_off_r])
                    
                else:
                    index_f = regex.search((primer_f4), sequence)
                    index_r = regex.search((primer_r4), sequence)
                    if index_f != None:
                        cut_off_f = int(index_f.span()[1])
                    if index_r != None:
                        cut_off_r = int(index_r.span()[0]+1)
                    if index_f or index_r != None:
                        tufA_list.append(record [cut_off_f:cut_off_r])
                    
                    else:
                        unmatched.append(record)
      
             
        
            
#Iterate over fastq file
print ("separating genes... It can take a while...")

original_reads = SeqIO.parse(input_file, "fastq")
do_it = primer_finder(original_reads,p16S_f,p16S_r,p18S_f,p18S_r,pUPA_f,pUPA_r,ptufA_f,ptufA_r)


count_um = SeqIO.write(unmatched, output_file_umatched, "fastq")
print ""
print "%i sequences did not match the primers and were stored in the %s file." %(count_um, output_file_umatched)
print ""

count_16S = SeqIO.write(s16_list, output_file_16S, "fastq")
print"",
print "Saved %i reads in the %s file." %(count_16S, output_file_16S)

count_18S = SeqIO.write(s18_list, output_file_18S, "fastq")
print""
print "Saved %i reads in the %s file." %(count_18S, output_file_18S)

count_UPA = SeqIO.write(UPA_list, output_file_UPA, "fastq")
print""
print "Saved %i reads in the %s file." %(count_UPA, output_file_UPA)

count_tufA = SeqIO.write(tufA_list, output_file_tufA, "fastq")
print""
print "Saved %i reads in the %s file." %(count_tufA, output_file_tufA)
print""
print "Done!"


# -*- coding: utf-8 -*-
"""
Script to separate genes from MiSeq run fastQ files (NEXTERA style)
AFTER the reads have been trimmed for the tail and adapter (see trim_adaptors.py)
It will produce one file for each gene (4 here - 16S, 18S, UPA, tufA)

Created on Mon Jul 28 17:27:33 2014
@author: VanessaRM
"""
""" USAGE:
give the fil name and tell if it is forward or reverse
ex:
python speratae_genes.py my_fastq_file.fastq forward
"""


from Bio import SeqIO
import regex
import sys

# Help text
if len(sys.argv) == 1:
    print ""
    print "Script to separate genes from a fastQ file"
    print ""
    print "Usage: supply the file name"
    print "ex: python speratae_genes.py my_fastq_file.fastq"
    print ""
    print "If you want to run it for multiple files, use the shell:"
    print "for file in *_R1_001.fastq; do python separate_genes.py $file; done >> screen.out 2>> screen.err &"
    print ""
    sys.exit()



# input files and arguments:
input_file = str(sys.argv[1])
output_file_16S = input_file + "_16S_.fastq"
output_file_18S = input_file + "_18.fastq"
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

unmatched_1=[] #seqs that didn't match with primers 

# Forward Primers (\w = ambiguous positions, only 12 "internal" bases used), allowing for 1 mismatch.
p16S_f = '(GTGCCAGC\wGCCGCGGTAA){e<=1}' #F515
p18S_f = '(GGTGGTGCATGGCCGTTCTTAGTT){e<=1}' #NF1
pUPA_f = '(GGACAGAAAGACCCTATGAA){e<=1}' #p23SrV_f1
ptufA_f = '(AC\wGG\wCG\wGG\wAC\wGT){e<=1}' #Oq_tuf2F

#Reverse primers(their reverse complement):
p16S_r = '(ATTAGA\wACCC\w\wGTAGTCC){e<=1}' #R806 rc
p18S_r = '(ATTACGTCCCTGCCCTTTGTA){e<=1}' # B18Sr2B rc
pUPA_r = '(CTCTAGGGATAACAGGCTGA){e<=1}' #p23SrV_r1 rc
ptufA_r = '(GCG\wTT\wGC\wATTCG\wGAAGG){e<=1}' #tufAR rc


# Define function primer_finder
def primer_finder (primer_f):
    global finder
    finder_fw = regex.search((primer_f), sequence)
    if finder_fw != None:
        finder = True
    else:
        finder = False
    return (finder)
    
    
#Iterate over fastq file
print ("separating genes... It can take more than an hour...")

# Run for forward primers:
for seq_record in SeqIO.parse(input_file , "fastq"):
    sequence = str(seq_record.seq)
    #    all_seqs.append(seq_record)
    
    primer_finder(p16S_f)
    if finder == True:
        s16_list.append (seq_record)       
    else:
            
            primer_finder(p18S_f)
            if finder == True:
                s18_list.append (seq_record)  
            else:
                    
                primer_finder(pUPA_f)
                if finder == True:
                    UPA_list.append (seq_record)
                else:
                                
                    primer_finder(ptufA_f)
                    if finder == True:
                        tufA_list.append (seq_record)                  
                    else:
                                 
# Run for reverse primers                                
                        primer_finder(p16S_r)
                        if finder == True:
                            s16_list.append (seq_record)       
                        else:
                            
                                primer_finder(p18S_r)
                                if finder == True:
                                    s18_list.append (seq_record)  
                                else:
                       
                                    primer_finder(pUPA_r)
                                    if finder == True:
                                        UPA_list.append (seq_record)
                                    else:
                               
                                        primer_finder(ptufA_r)
                                        if finder == True:
                                            tufA_list.append (seq_record)                  
                                        else:
                                            unmatched_1.append(seq_record)


#Save files
SeqIO.write(s16_list, output_file_16S, "fastq")
SeqIO.write(s18_list, output_file_18S, "fastq")
SeqIO.write(UPA_list, output_file_UPA, "fastq")
SeqIO.write(tufA_list, output_file_tufA, "fastq")
SeqIO.write(unmatched_1, output_file_umatched, "fastq")


#print results
print ("")
print ("%i sequences did not match the primers and were stored in the unmatched file" %(len(unmatched_1)))
print ("")
print ("%i sequences were stored in 16S file" %(len(s16_list)))
print ("%i sequences were stored in 18S file" %(len(s18_list)))
print ("%i sequences were stored in UPA file" %(len(UPA_list)))
print ("%i sequences were stored in TufA file" %(len(tufA_list)))
print ("")
print ("Done!")



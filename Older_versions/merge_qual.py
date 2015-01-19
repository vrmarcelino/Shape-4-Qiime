'''
Merge_quality_scores
add a random quality to the barcode read
@author: Vanessa
'''
barcode_qual = "25 25 25 25 25 25 25 25 25 25 25 25 "

from Bio import SeqIO
import re

adapted_qual =[]

for seq_record in SeqIO.parse("tufA_001_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)


for seq_record in SeqIO.parse("tufA_002_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)


for seq_record in SeqIO.parse("tufA_003_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)


for seq_record in SeqIO.parse("tufA_004_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)

for seq_record in SeqIO.parse("tufA_005_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)


for seq_record in SeqIO.parse("tufA_006_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)


for seq_record in SeqIO.parse("tufA_007_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)

for seq_record in SeqIO.parse("tufA_008_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)

for seq_record in SeqIO.parse("tufA_009_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)

for seq_record in SeqIO.parse("tufA_010_trim.qual", "qual"):
    original_qual = str(seq_record.letter_annotations["phred_quality"])
    new_qual = (barcode_qual + original_qual)
    new_qual = re.sub(r"\[|\]|\,","",new_qual) # delete [ and commas
    name = ">" + (seq_record.description)
    
    adapted_qual.append (name)
    adapted_qual.append (new_qual)

#save
savefile = open("tufA_concat.qual", "w")
for lines in adapted_qual:
    savefile.write("%s\n" % lines)
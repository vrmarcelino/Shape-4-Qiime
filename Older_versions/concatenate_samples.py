'''
Concatenate the fasta files (different samples/same gene fragment)
    for downstream analysis in QIIME 
'''
from Bio import SeqIO

# Define Barcodes (Earth microbiome project - Golay barcodes -806rcbc0 - 806rcbc9)
barcode_001 = "TCCCTTGTCTCC"
barcode_002 = "ACGAGACTGATT"
barcode_003 = "GCTGTACGGATT"
barcode_004 = "ATCACCAGGTGT"
barcode_005 = "TGGTCAACGATA"
barcode_006 = "ATCGCACAGTAA"
barcode_007 = "GTCGTGTAGCCT"
barcode_008 = "AGCGGAGGTTAG"
barcode_009 = "ATCCTTTGGTTC"
barcode_010 = "TACAGCGCATAC"

#Store the files
all_records = []

print ("parsing file 001...")
for seq_record in SeqIO.parse("fasta_con_1.fas", "fasta"):
    seq_record.seq = (barcode_001 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 002...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_002 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 003...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_003 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 004...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_004 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 005...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_005 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 006...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_006 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 007...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_007 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 008...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_008 + seq_record.seq)
    all_records.append (seq_record)

print ("parsing file 009...")
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_009 + seq_record.seq)
    all_records.append (seq_record)
    
print ("parsing last file...")   
for seq_record in SeqIO.parse("fasta_con_2.fas", "fasta"):
    seq_record.seq = (barcode_010 + seq_record.seq)
    all_records.append (seq_record)

SeqIO.write(all_records, "all_records.fna", "fasta")
print ("Done!")
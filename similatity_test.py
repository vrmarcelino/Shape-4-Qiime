# -*- coding: utf-8 -*-
"""
Check similarity among same species/genera sequences
USAGE: similarity_test.py -t ref_alignment.fasta -t species

Warning: this script simply computes the fraction of non-variant sites. 
It gives a rough idea of the threshold that you want to use to define OTUs, but ne aware that it will vary with 
the number of seqs each species has in the reference alignment. 
Here it is set to consider only taxa with more than 10 representatives.

@author: VanessaRM
Created on Thu Mar 12 15:12:41 2015

"""
#load packages:
from Bio import SeqIO
from Bio import AlignIO
from argparse import ArgumentParser
import os
import numpy as np

#Input data
parser = ArgumentParser()
parser.add_argument('-i', '--alignment', help='The path to the reference alignment in fasta format', required=True)
parser.add_argument('-t', '--tax_level', help='The taxonomic level you want to analyse - species or genera, DEFAULT: species', required=False)

args = parser.parse_args()
seqs = AlignIO.read(args.alignment, "fasta")
wish_tax_level = str(args.tax_level) #species or genus, species is default

seqs_by_group_dict = dict()
simil_frac = 0.0

# function to find the species name in fasta file and store seqs in different dictionaries
def tax_level_function (sequence_rec, wish_tax_level="species"):
    whole_seq_identification = sequence_rec.id
    name_parts = whole_seq_identification.split("_")
    genus_name = name_parts[0]
    species_name = name_parts[1]
    sequence = sequence_rec

#sort out tax level to group the seqs     
    if wish_tax_level == "species":
        group_factor = species_name
    elif wish_tax_level == "genera":
        group_factor = genus_name
    else:
        print "Define threshold by species or genus?"
      
    
# group the seqs in a dictionary  
    try:
        #append the new seq to the existing array for that species
        seqs_by_group_dict[group_factor].append (sequence)
        
    except KeyError:
        #creates a new array in this slot
        seqs_by_group_dict[group_factor] = [sequence]
        
    return seqs_by_group_dict


# Function to calculate the similarity fraction
# note that the similarity will be lower the more sequences you have in the alignment    
def sim_frac_calculation (alignemnt):
    global simil_frac
    count_zero_variation = 0
    align_lengh = alignemnt.get_alignment_length()
    for r in range (0, align_lengh):
        snps_with_gaps = ''.join(set(alignemnt[:,r]))
        snps = snps_with_gaps.replace("-", "") # note that sites with a gap and a base will count as 1 snps!
    
        #fraction of sites where snps = 1:
        if len(snps) > 1: # = 2 or more
            count_zero_variation += 1

        simil_frac = 1 - (float(count_zero_variation) / float(align_lengh))
    return simil_frac
    

# Run the function to divide the seqs into groups by species or genus
# Will return seqs_by_group_dict
for s in seqs:
    tax_level_function (s, wish_tax_level)
    

store_similarities = []
for group_id in seqs_by_group_dict:

#exclude alignemnts with less than 10 sequences
    if len (seqs_by_group_dict[group_id]) > 9:
                    
        # convert to alignment again            
        species_aln = SeqIO.write(seqs_by_group_dict[group_id],group_id,"fasta")
        x_aln = AlignIO.read(group_id, "fasta")
        sim_frac_calculation(x_aln)
        store_similarities.append (simil_frac)
        print group_id
        print simil_frac
        print ""
        
        os.remove(group_id)
        
mean_similarity = np.mean(store_similarities)
print "The mean similarity among these %s is %s" %(wish_tax_level, round(mean_similarity,2))




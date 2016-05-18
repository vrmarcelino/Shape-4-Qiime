# -*- coding: utf-8 -*-
# delete seqs from align that are not present in a tree

os.chdir('../15_phylobayes')

from Bio import SeqIO
from Bio import Phylo


tree_file_path = "new_tree.nwk" # substitute by ArgumentParser
tree = Phylo.read(tree_file_path, "newick")

new_alignment = []
counter = 0
for seq_record in SeqIO.parse("final_aln.fas" , "fasta"):
    if seq_record.id in str(tree.get_terminals()):
        new_alignment.append(seq_record)
    else:
        counter += 1

SeqIO.write(new_alignment, "new_align.fas", "fasta")

print "%i taxa were excluded" %(counter)

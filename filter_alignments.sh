#!/bin/bash

echo "running trimAl"
### Run trimAl to filter spurious sequences from mafft alignments.
trimal -in 05_alignments/01_intermediates/UPA.nt_aligned.fas.reduced_pfiltered.fasta -out 05_alignments/01_intermediates/UPA_trimal.fasta -resoverlap 0.80 -seqoverlap 80 >> screen_UPA.out 2>> screen_UPA.err
trimal -in 05_alignments/01_intermediates/tufA.nt_aligned.fas.reduced_pfiltered.fasta -out 05_alignments/01_intermediates/tufA_trimal.fasta -resoverlap 0.80 -seqoverlap 80 >> screen_tufA.out 2>> screen_tufA.err

echo "running get_filtered_otus.py"
### Run get_filtered_otus.py to obtain a fasta file with sequences that didn't align properly
get_filtered.otus.py -oa 05_alignments/01_intermediates/UPA.nt_aligned.fas.reduced_pfiltered.fasta -fa 05_alignments/01_intermediates/UPA_trimal.fasta  -o 05_alignments/01_intermediates/UPA.otus_failures.fasta >> screen_UPA.out 2>> screen_UPA.err
get_filtered.otus.py -oa 05_alignments/01_intermediates/tufA.nt_aligned.fas.reduced_pfiltered.fasta -fa 05_alignments/01_intermediates/tufA_trimal.fasta  -o 05_alignments/01_intermediates/tufA.otus_failures.fasta >> screen_tufA.out 2>> screen_tufA.err

echo "done"
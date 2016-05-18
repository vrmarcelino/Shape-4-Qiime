#!/bin/bash

### Step 10
### OTU filtering - continuation
### It varies according to the study

echo "selecting study-specififc OTUs"
# Filter samples with zero counts, controls and samples that not belong to the diversity paper:
filter_samples_from_otu_table.py -i 07_OTU_table/16S/OTU_table_3.biom -o 07_OTU_table/16S/OTU_table_4.biom -m Mapping/map.txt -s "tufDiversity:keep" -n 1
filter_samples_from_otu_table.py -i 07_OTU_table/18S/OTU_table_3.biom -o 07_OTU_table/18S/OTU_table_4.biom -m Mapping/map.txt -s "tufDiversity:keep" -n 1 
filter_samples_from_otu_table.py -i 07_OTU_table/UPA/OTU_table_3.biom -o 07_OTU_table/UPA/OTU_table_4.biom -m Mapping/map.txt -s "tufDiversity:keep" -n 1
filter_samples_from_otu_table.py -i 07_OTU_table/tufA/OTU_table_3.biom -o 07_OTU_table/tufA/OTU_table_4.biom -m Mapping/map.txt -s "tufDiversity:keep" -n 1

# Print third statistics
biom summarize_table -i 07_OTU_table/16S/OTU_table_4.biom -o 07_OTU_table/16S/OTUs_stats3.txt
biom summarize_table -i 07_OTU_table/18S/OTU_table_4.biom -o 07_OTU_table/18S/OTUs_stats3.txt
biom summarize_table -i 07_OTU_table/UPA/OTU_table_4.biom -o 07_OTU_table/UPA/OTUs_stats3.txt
biom summarize_table -i 07_OTU_table/tufA/OTU_table_4.biom -o 07_OTU_table/tufA/OTUs_stats3.txt

echo "Producing community composition statistics... it takes a while..."
# Taxonomic stats
mkdir 08_Tax_Composition
summarize_taxa_through_plots.py -i 07_OTU_table/16S/OTU_table_4.biom -o 08_Tax_Composition/16S -m Mapping/map.txt
summarize_taxa_through_plots.py -i 07_OTU_table/18S/OTU_table_4.biom -o 08_Tax_Composition/18S -m Mapping/map.txt
summarize_taxa_through_plots.py -i 07_OTU_table/UPA/OTU_table_4.biom -o 08_Tax_Composition/UPA -m Mapping/map.txt
summarize_taxa_through_plots.py -i 07_OTU_table/tufA/OTU_table_4.biom -o 08_Tax_Composition/tufA -m Mapping/map.txt

echo "Getting eukaryitic OTUs"
# Filter by taxonomy: only d_Eukaryota for phylonegy (search - Chloroplast for 16S)
filter_taxa_from_otu_table.py -i 07_OTU_table/16S/OTU_table_4.biom -p c__Chloroplast -o 07_OTU_table/16S/OTU_table_5.biom
filter_taxa_from_otu_table.py -i 07_OTU_table/18S/OTU_table_4.biom -p Eukaryota -o 07_OTU_table/18S/OTU_table_5.biom
filter_taxa_from_otu_table.py -i 07_OTU_table/UPA/OTU_table_4.biom -p d_Eukaryota -o 07_OTU_table/UPA/OTU_table_5.biom
filter_taxa_from_otu_table.py -i 07_OTU_table/tufA/OTU_table_4.biom -p d_Eukaryota -o 07_OTU_table/tufA/OTU_table_5.biom

# Filter samples again to exclude samples that didn't have any Eukaryote:
filter_samples_from_otu_table.py -i 07_OTU_table/16S/OTU_table_5.biom -o 07_OTU_table/16S/OTU_table_6.biom -n 1
filter_samples_from_otu_table.py -i 07_OTU_table/18S/OTU_table_5.biom -o 07_OTU_table/18S/OTU_table_6.biom -n 1
filter_samples_from_otu_table.py -i 07_OTU_table/UPA/OTU_table_5.biom -o 07_OTU_table/UPA/OTU_table_6.biom -n 1
filter_samples_from_otu_table.py -i 07_OTU_table/tufA/OTU_table_5.biom -o 07_OTU_table/tufA/OTU_table_6.biom -n 1

# Print fourth statistics
biom summarize_table -i 07_OTU_table/16S/OTU_table_6.biom -o 07_OTU_table/16S/OTUs_stats4.txt
biom summarize_table -i 07_OTU_table/18S/OTU_table_6.biom -o 07_OTU_table/18S/OTUs_stats4.txt
biom summarize_table -i 07_OTU_table/UPA/OTU_table_6.biom -o 07_OTU_table/UPA/OTUs_stats4.txt
biom summarize_table -i 07_OTU_table/tufA/OTU_table_6.biom -o 07_OTU_table/tufA/OTUs_stats4.txt

### Filter reference sequences:
echo "Filtering sequences..."
mkdir 09_Euk_Otus
filter_fasta.py -f 04_blastfilter/02_final_results/16S.otus.fasta -b 07_OTU_table/16S/OTU_table_6.biom -o 09_Euk_Otus/16S.euk_otus.fas
filter_fasta.py -f 04_blastfilter/02_final_results/18S.otus.fasta -b 07_OTU_table/18S/OTU_table_6.biom -o 09_Euk_Otus/18S.euk_otus.fas
filter_fasta.py -f 04_blastfilter/02_final_results/UPA.otus.fasta -b 07_OTU_table/UPA/OTU_table_6.biom -o 09_Euk_Otus/UPA.euk_otus.fas
filter_fasta.py -f 04_blastfilter/02_final_results/tufA.otus.fasta -b 07_OTU_table/tufA/OTU_table_6.biom -o 09_Euk_Otus/tufA.euk_otus.fas

echo "Filtering aligned sequences"
filter_fasta.py -f 05_alignments/01_intermediates/16S_aligned/16S.otus_aligned_pfiltered.fasta -b 07_OTU_table/16S/OTU_table_6.biom -o 09_Euk_Otus/16S.euk_otus_aln.fas
filter_fasta.py -f 05_alignments/01_intermediates/18S_aligned/18S.otus_aligned_pfiltered.fasta -b 07_OTU_table/18S/OTU_table_6.biom -o 09_Euk_Otus/18S.euk_otus_aln.fas
filter_fasta.py -f 05_alignments/01_intermediates/UPA_trimal.fasta -b 07_OTU_table/UPA/OTU_table_6.biom -o 09_Euk_Otus/UPA.euk_otus_aln.fas
filter_fasta.py -f 05_alignments/01_intermediates/tufA_trimal.fasta -b 07_OTU_table/tufA/OTU_table_6.biom -o 09_Euk_Otus/tufA.euk_otus_aln.fas


echo "Producing community composition statistics with Eukyotic OTUs only... "
# Taxonomic stats
mkdir 10_Euk_Tax_Composition
summarize_taxa_through_plots.py -i 07_OTU_table/16S/OTU_table_6.biom -o 10_Euk_Tax_Composition/16S -m Mapping/map.txt
summarize_taxa_through_plots.py -i 07_OTU_table/18S/OTU_table_6.biom -o 10_Euk_Tax_Composition/18S -m Mapping/map.txt
summarize_taxa_through_plots.py -i 07_OTU_table/UPA/OTU_table_6.biom -o 10_Euk_Tax_Composition/UPA -m Mapping/map.txt
summarize_taxa_through_plots.py -i 07_OTU_table/tufA/OTU_table_6.biom -o 10_Euk_Tax_Composition/tufA -m Mapping/map.txt

echo "Done"



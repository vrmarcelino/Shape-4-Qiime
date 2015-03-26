#!/bin/bash

# Pipeline for analysing environmental samples.
# Part 3 - Clean dataset from spurious reads, make and clean OTU table
# Part of this pipeline only works with Qiime 1.8
# Use make_otu_heatmap.py instead of make_otu_heatmap_html.py if using Qiime 1.9
# Requires patchines_ids.txt or other txt file with samples IDs to keep in the following analyses
# 11 - Feb - 2015

echo "Launching....."

#only need if running macqiime
source /macqiime/configs/bash_profile.txt


### Make OTU_table.biom and filter OTUs with less than 5 reads (in total).
# 16S - already has an OTU_table.biom, but it need to be rebuilt to work with Qiime 1.8
# 18S - needs to include taxonomy
# Filter pynast alignmnet failures

echo "Making OTU tables"

mkdir 10_OTU_filtering
mkdir 10_OTU_filtering/16S
mkdir 10_OTU_filtering/18S
mkdir 10_OTU_filtering/UPA
mkdir 10_OTU_filtering/tufA

make_otu_table.py -i 08_OTUs/16S/final_otu_map_mc5.txt -t 08_OTUs/16S/uclust_assigned_taxonomy/rep_set_tax_assignments.txt -o 10_OTU_filtering/16S/otu_table_mc5_w_tax.biom >> screen_16S.out 2>> screen_16S.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/16S/otu_table_mc5_w_tax.biom -o 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures.biom --min_count 5 -e 09_Alignment/16S/rep_set_failures.fasta >> screen_16S.out 2>> screen_16S.err

make_otu_table.py -i 08_OTUs/18S/final_otu_map_mc5.txt -t 08_OTUs/18S/uclust_assigned_taxonomy/rep_set_tax_assignments.txt -o 10_OTU_filtering/18S/otu_table_mc5_w_tax.biom >> screen_18S.out 2>> screen_18S.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/18S/otu_table_mc5_w_tax.biom -o 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures.biom --min_count 5 -e 09_Alignment/18S/rep_set_failures.fasta >> screen_18S.out 2>> screen_18S.err

make_otu_table.py -i 08_OTUs/UPA/seqs_no_chimeras_otus.txt -t 08_OTUs/UPA/rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o 10_OTU_filtering/UPA/OTU_table.biom >> screen_UPA.out 2>> screen_UPA.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/UPA/OTU_table.biom -o 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures.biom --min_count 5 -e 09_Alignment/UPA/rep_set_failures.fasta >> screen_UPA.out 2>> screen_UPA.err

make_otu_table.py -i 08_OTUs/tufA/seqs_no_chimeras_otus.txt -t 08_OTUs/tufA/rdp_assigned_taxonomy/rep_set_tax_assignments.txt -o 10_OTU_filtering/tufA/OTU_table.biom >> screen_tufA.out 2>> screen_tufA.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/tufA/OTU_table.biom -o 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures.biom --min_count 5 -e 09_Alignment/tufA/rep_set_failures.fasta >> screen_tufA.out 2>> screen_tufA.err


### Filter OTU table by min. count of reads per OTU per sample. = 10 
# Requires filter_observations_by_sample.py (https://gist.github.com/adamrp/7591573). Adam did a script to work with Qiime 1.9 as well.
echo "Filtering by sample"
filter_observations_by_sample.py -i 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures.biom -n 10 -o 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom >> screen_16S.out 2>> screen_16S.err
filter_observations_by_sample.py -i 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures.biom -n 10 -o 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom >> screen_18S.out 2>> screen_18S.err
filter_observations_by_sample.py -i 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures.biom -n 10 -o 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom >> screen_UPA.out 2>> screen_UPA.err
filter_observations_by_sample.py -i 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures.biom -n 10 -o 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom >> screen_tufA.out 2>> screen_tufA.err


echo "Filtering control samples"
### Filter OTUs found in control samples

#create a mapping file for controls (see mapping file)
filter_samples_from_otu_table.py -i 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/16S/control_samples.biom -m 05_Mapping/16S/Map_16S.txt -s "Host:control" >> screen_16S.out 2>> screen_16S.err
filter_samples_from_otu_table.py -i 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/18S/control_samples.biom -m 05_Mapping/18S/Map_18S.txt -s "Host:control" >> screen_18S.out 2>> screen_18S.err
filter_samples_from_otu_table.py -i 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/UPA/control_samples.biom -m 05_Mapping/UPA/Map_UPA.txt -s "Host:control" >> screen_UPA.out 2>> screen_UPA.err
filter_samples_from_otu_table.py -i 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/tufA/control_samples.biom -m 05_Mapping/tufA/Map_tufA.txt -s "Host:control" >> screen_tufA.out 2>> screen_tufA.err

# filter OTUs with 0 counts
filter_otus_from_otu_table.py -i 10_OTU_filtering/16S/control_samples.biom -o 10_OTU_filtering/16S/control_samples_filtered.biom -n 1 >> screen_16S.out 2>> screen_16S.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/18S/control_samples.biom -o 10_OTU_filtering/18S/control_samples_filtered.biom -n 1 >> screen_18S.out 2>> screen_18S.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/UPA/control_samples.biom -o 10_OTU_filtering/UPA/control_samples_filtered.biom -n 1 >> screen_UPA.out 2>> screen_UPA.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/tufA/control_samples.biom -o 10_OTU_filtering/tufA/control_samples_filtered.biom -n 1 >> screen_tufA.out 2>> screen_tufA.err

#convert to txt file: 
biom convert -b -i 10_OTU_filtering/16S/control_samples_filtered.biom -o 10_OTU_filtering/16S/otus_to_remove.txt >> screen_16S.out 2>> screen_16S.err
biom convert -b -i 10_OTU_filtering/18S/control_samples_filtered.biom -o 10_OTU_filtering/18S/otus_to_remove.txt >> screen_18S.out 2>> screen_18S.err
biom convert -b -i 10_OTU_filtering/UPA/control_samples_filtered.biom -o 10_OTU_filtering/UPA/otus_to_remove.txt >> screen_UPA.out 2>> screen_UPA.err
biom convert -b -i 10_OTU_filtering/tufA/control_samples_filtered.biom -o 10_OTU_filtering/tufA/otus_to_remove.txt >> screen_tufA.out 2>> screen_tufA.err

# Remove only OTUs representing less than 2% of the reads:
# Requires all_otus.txt + otus_to_remove.txt and filter_otus_to_exclude.py
# get all otus.txt
biom convert -b -i 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/16S/all_otus.txt >> screen_16S.out 2>> screen_16S.err
biom convert -b -i 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/18S/all_otus.txt >> screen_18S.out 2>> screen_18S.err
biom convert -b -i 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/UPA/all_otus.txt >> screen_UPA.out 2>> screen_UPA.err
biom convert -b -i 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/tufA/all_otus.txt >> screen_tufA.out 2>> screen_tufA.err
# Filter real contamninats by 2% 
filter_otus_to_exclude.py -e 10_OTU_filtering/16S/otus_to_remove.txt -a 10_OTU_filtering/16S/all_otus.txt -o 10_OTU_filtering/16S/otus_to_really_exclude.txt -t 2 >> screen_16S.out 2>> screen_16S.err
filter_otus_to_exclude.py -e 10_OTU_filtering/18S/otus_to_remove.txt -a 10_OTU_filtering/18S/all_otus.txt -o 10_OTU_filtering/18S/otus_to_really_exclude.txt -t 2 >> screen_18S.out 2>> screen_18S.err
filter_otus_to_exclude.py -e 10_OTU_filtering/UPA/otus_to_remove.txt -a 10_OTU_filtering/UPA/all_otus.txt -o 10_OTU_filtering/UPA/otus_to_really_exclude.txt -t 2 >> screen_UPA.out 2>> screen_UPA.err
filter_otus_to_exclude.py -e 10_OTU_filtering/tufA/otus_to_remove.txt -a 10_OTU_filtering/tufA/all_otus.txt -o 10_OTU_filtering/tufA/otus_to_really_exclude.txt -t 2 >> screen_tufA.out 2>> screen_tufA.err
# Filter otus based on the new otus_to_really_remove table
filter_otus_from_otu_table.py -i 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -e 10_OTU_filtering/16S/otus_to_really_exclude.txt >> screen_16S.out 2>> screen_16S.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -e 10_OTU_filtering/18S/otus_to_really_exclude.txt >> screen_18S.out 2>> screen_18S.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -e 10_OTU_filtering/UPA/otus_to_really_exclude.txt >> screen_UPA.out 2>> screen_UPA.err
filter_otus_from_otu_table.py -i 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts.biom -o 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -e 10_OTU_filtering/tufA/otus_to_really_exclude.txt >> screen_tufA.out 2>> screen_tufA.err

# Filter samples with zero counts, controls and samples that not belong to the patchiness experiment:
#filter_samples_from_otu_table.py -i 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/16S/otu_table_final.biom  --sample_id_fp acidification_ids.txt -n 1 >> screen_16S.out 2>> screen_16S.err
#filter_samples_from_otu_table.py -i 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/18S/otu_table_final.biom --sample_id_fp acidification_ids.txt -n 1 >> screen_18S.out 2>> screen_18S.err
#filter_samples_from_otu_table.py -i 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/UPA/otu_table_final.biom --sample_id_fp acidification_ids.txt -n 1 >> screen_UPA.out 2>> screen_UPA.err
#filter_samples_from_otu_table.py -i 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/tufA/otu_table_semi_final.biom --sample_id_fp acidification_ids.txt -n 1 >> screen_tufA.out 2>> screen_tufA.err

# filter only samples with zero counts
filter_samples_from_otu_table.py -i 10_OTU_filtering/16S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/16S/otu_table_final.biom -n 1 >> screen_16S.out 2>> screen_16S.err
filter_samples_from_otu_table.py -i 10_OTU_filtering/18S/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/18S/otu_table_final.biom -n 1 >> screen_18S.out 2>> screen_18S.err
filter_samples_from_otu_table.py -i 10_OTU_filtering/UPA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/UPA/otu_table_final.biom -n 1 >> screen_UPA.out 2>> screen_UPA.err
filter_samples_from_otu_table.py -i 10_OTU_filtering/tufA/otu_table_mc5_w_tax_no_pynast_failures_no_ghosts_no_contam.biom -o 10_OTU_filtering/tufA/otu_table_semi_final.biom -n 1 >> screen_tufA.out 2>> screen_tufA.err


# Filter by taxonomy: only Ulvophyceae for tufA dataset
filter_taxa_from_otu_table.py -i 10_OTU_filtering/tufA/otu_table_semi_final.biom -p c_Ulvophyceae -o 10_OTU_filtering/tufA/otu_table_final.biom >> screen_tufA.out 2>> screen_tufA.err



### Check with OTU stats and heatmap:
echo "Making stats and heatmaps"

biom summarize_table -i 10_OTU_filtering/16S/otu_table_final.biom -o 10_OTU_filtering/16S/OTUs_stats.txt >> screen_16S.out 2>> screen_16S.err
make_otu_heatmap_html.py -i 10_OTU_filtering/16S/otu_table_final.biom -o 10_OTU_filtering/16S/OTU_Heatmap >> screen_16S.out 2>> screen_16S.err

biom summarize_table -i 10_OTU_filtering/18S/otu_table_final.biom -o 10_OTU_filtering/18S/OTUs_stats.txt >> screen_18S.out 2>> screen_18S.err
make_otu_heatmap_html.py -i 10_OTU_filtering/18S/otu_table_final.biom -o 10_OTU_filtering/18S/OTU_Heatmap >> screen_18S.out 2>> screen_18S.err

biom summarize_table -i 10_OTU_filtering/UPA/otu_table_final.biom -o 10_OTU_filtering/UPA/OTUs_stats.txt >> screen_UPA.out 2>> screen_UPA.err
make_otu_heatmap_html.py -i 10_OTU_filtering/UPA/otu_table_final.biom -o 10_OTU_filtering/UPA/OTU_Heatmap >> screen_UPA.out 2>> screen_UPA.err

biom summarize_table -i 10_OTU_filtering/tufA/otu_table_final.biom -o 10_OTU_filtering/tufA/OTUs_stats.txt >> screen_tufA.out 2>> screen_tufA.err
make_otu_heatmap_html.py -i 10_OTU_filtering/tufA/otu_table_final.biom -o 10_OTU_filtering/tufA/OTU_Heatmap >> screen_tufA.out 2>> screen_tufA.err

### Filter reference sequences:
echo "Filtering Reference sequences..."
filter_fasta.py -f 08_OTUs/16S/rep_set.fna -b 10_OTU_filtering/16S/otu_table_final.biom -o 08_OTUs/16S/rep_set_final.fna >> screen_16S.out 2>> screen_16S.err
filter_fasta.py -f 08_OTUs/18S/rep_set.fna -b 10_OTU_filtering/18S/otu_table_final.biom -o 08_OTUs/18S/rep_set_final.fna >> screen_18S.out 2>> screen_18S.err
filter_fasta.py -f 08_OTUs/UPA/rep_set.fna -b 10_OTU_filtering/UPA/otu_table_final.biom -o 08_OTUs/UPA/rep_set_final.fna >> screen_UPA.out 2>> screen_UPA.err
filter_fasta.py -f 08_OTUs/tufA/rep_set.fna -b 10_OTU_filtering/tufA/otu_table_final.biom -o 08_OTUs/tufA/rep_set_final.fna >> screen_tufA.out 2>> screen_tufA.err

echo "Filtering aligned reference sequences..."
filter_fasta.py -f 09_Alignment/16S/rep_set_aligned_pfiltered.fasta -b 10_OTU_filtering/16S/otu_table_final.biom -o 09_Alignment/16S/rep_set_aligned_biomfiltered.fasta >> screen_16S.out 2>> screen_16S.err
filter_fasta.py -f 09_Alignment/18S/rep_set_aligned_pfiltered.fasta -b 10_OTU_filtering/18S/otu_table_final.biom -o 09_Alignment/18S/rep_set_aligned_biomfiltered.fasta >> screen_18S.out 2>> screen_168.err
filter_fasta.py -f 09_Alignment/UPA/rep_set_aligned_pfiltered.fasta -b 10_OTU_filtering/UPA/otu_table_final.biom -o 09_Alignment/UPA/rep_set_aligned_biomfiltered.fasta >> screen_UPA.out 2>> screen_UPA.err
filter_fasta.py -f 09_Alignment/tufA/rep_set_aligned_pfiltered.fasta -b 10_OTU_filtering/tufA/otu_table_final.biom -o 09_Alignment/tufA/rep_set_aligned_biomfiltered.fasta >> screen_tufA.out 2>> screen_tufA.err


echo "Done"
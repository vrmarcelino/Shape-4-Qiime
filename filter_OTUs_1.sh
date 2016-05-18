#!/bin/bash

### Step 9
### OTU filtering

mkdir 07_OTU_table
mkdir 07_OTU_table/16S
mkdir 07_OTU_table/18S
mkdir 07_OTU_table/UPA
mkdir 07_OTU_table/tufA


#Convert UC (otu map from uparse) to otu-table.txt (qiime compatible format) <<< Brazilian Microbiome Project SCRIPT: bmp-map2qiime.py>>>
echo "converting OTU table to Qiime format"
bmp-map2qiime.py 04_blastfilter/02_final_results/16S.map.uc > 07_OTU_table/16S/otu_table.txt
bmp-map2qiime.py 04_blastfilter/02_final_results/18S.map.uc > 07_OTU_table/18S/otu_table.txt
bmp-map2qiime.py 04_blastfilter/02_final_results/UPA.map.uc > 07_OTU_table/UPA/otu_table.txt
bmp-map2qiime.py 04_blastfilter/02_final_results/tufA.map.uc > 07_OTU_table/tufA/otu_table.txt


### Make otu_table
echo "Filtering otus that didn't aligned, mapping taxonomy and producing first stats"
# Make biom table and include taxonomic assignments
make_otu_table.py -i 07_OTU_table/16S/otu_table.txt -t 06_Taxonomy/16S/16S.otus_tax_assignments.txt -o 07_OTU_table/16S/OTU_table_0.biom
make_otu_table.py -i 07_OTU_table/18S/otu_table.txt -t 06_Taxonomy/18S/18S.otus_tax_assignments.txt -o 07_OTU_table/18S/OTU_table_0.biom
make_otu_table.py -i 07_OTU_table/UPA/otu_table.txt -t 06_Taxonomy/UPA/UPA.otus_tax_assignments.txt -o 07_OTU_table/UPA/OTU_table_0.biom
make_otu_table.py -i 07_OTU_table/tufA/otu_table.txt -t 06_Taxonomy/tufA/tufA.otus_tax_assignments.txt -o 07_OTU_table/tufA/OTU_table_0.biom

# Filter OTUs that didn't align and otus with less than 5 reads.
filter_otus_from_otu_table.py -i 07_OTU_table/16S/OTU_table_0.biom -o 07_OTU_table/16S/OTU_table_1.biom --min_count 5 -e 05_alignments/01_intermediates/16S_aligned/16S.otus_failures.fasta 
filter_otus_from_otu_table.py -i 07_OTU_table/18S/OTU_table_0.biom -o 07_OTU_table/18S/OTU_table_1.biom --min_count 5 -e 05_alignments/01_intermediates/18S_aligned/18S.otus_failures.fasta 
filter_otus_from_otu_table.py -i 07_OTU_table/UPA/OTU_table_0.biom -o 07_OTU_table/UPA/OTU_table_1.biom --min_count 5 -e 05_alignments/01_intermediates/UPA.otus_failures.fasta 
filter_otus_from_otu_table.py -i 07_OTU_table/tufA/OTU_table_0.biom -o 07_OTU_table/tufA/OTU_table_1.biom --min_count 5 -e 05_alignments/01_intermediates/tufA.otus_failures.fasta

# Print first statistics
biom summarize_table -i 07_OTU_table/16S/OTU_table_1.biom -o 07_OTU_table/16S/OTUs_stats1.txt
biom summarize_table -i 07_OTU_table/18S/OTU_table_1.biom -o 07_OTU_table/18S/OTUs_stats1.txt
biom summarize_table -i 07_OTU_table/UPA/OTU_table_1.biom -o 07_OTU_table/UPA/OTUs_stats1.txt
biom summarize_table -i 07_OTU_table/tufA/OTU_table_1.biom -o 07_OTU_table/tufA/OTUs_stats1.txt


### Filter OTU table by min. count of reads per OTU per sample. = 2
# You may want to increase the threshold if doing beta diveristy analysis. For alpha diversity it doesn't matter.
# Requires filter_observations_by_sample.py (https://gist.github.com/adamrp/7591573). Adam did a script to work with Qiime 1.9 as well.
echo "Filtering by sample - min 2 reads per otu per sample "
filter_observations_by_sample.py -i 07_OTU_table/16S/OTU_table_1.biom -n 2 -o 07_OTU_table/16S/OTU_table_2.biom
filter_observations_by_sample.py -i 07_OTU_table/18S/OTU_table_1.biom -n 2 -o 07_OTU_table/18S/OTU_table_2.biom
filter_observations_by_sample.py -i 07_OTU_table/UPA/OTU_table_1.biom -n 2 -o 07_OTU_table/UPA/OTU_table_2.biom
filter_observations_by_sample.py -i 07_OTU_table/tufA/OTU_table_1.biom -n 2 -o 07_OTU_table/tufA/OTU_table_2.biom


### Filter OTUs found in control samples
echo "Filtering controls"
#create a mapping file for controls 
### Filter OTUs found in control samples

#create a mapping file for controls:
filter_samples_from_otu_table.py -i 07_OTU_table/16S/OTU_table_2.biom -o 07_OTU_table/16S/control_samples.biom -m Mapping/map.txt -s "Host:control"
filter_samples_from_otu_table.py -i 07_OTU_table/18S/OTU_table_2.biom -o 07_OTU_table/18S/control_samples.biom -m Mapping/map.txt -s "Host:control"
filter_samples_from_otu_table.py -i 07_OTU_table/UPA/OTU_table_2.biom -o 07_OTU_table/UPA/control_samples.biom -m Mapping/map.txt -s "Host:control"
filter_samples_from_otu_table.py -i 07_OTU_table/tufA/OTU_table_2.biom -o 07_OTU_table/tufA/control_samples.biom -m Mapping/map.txt -s "Host:control"

# filter OTUs with 0 counts
filter_otus_from_otu_table.py -i 07_OTU_table/16S/control_samples.biom -o 07_OTU_table/16S/control_samples_filtered.biom -n 1
filter_otus_from_otu_table.py -i 07_OTU_table/18S/control_samples.biom -o 07_OTU_table/18S/control_samples_filtered.biom -n 1
filter_otus_from_otu_table.py -i 07_OTU_table/UPA/control_samples.biom -o 07_OTU_table/UPA/control_samples_filtered.biom -n 1
filter_otus_from_otu_table.py -i 07_OTU_table/tufA/control_samples.biom -o 07_OTU_table/tufA/control_samples_filtered.biom -n 1

#convert to txt file: 
biom convert -i 07_OTU_table/16S/control_samples_filtered.biom -o 07_OTU_table/16S/otus_to_remove.txt --to-tsv --table-type="OTU table"
biom convert -i 07_OTU_table/18S/control_samples_filtered.biom -o 07_OTU_table/18S/otus_to_remove.txt --to-tsv --table-type="OTU table"
biom convert -i 07_OTU_table/UPA/control_samples_filtered.biom -o 07_OTU_table/UPA/otus_to_remove.txt --to-tsv --table-type="OTU table"
biom convert -i 07_OTU_table/tufA/control_samples_filtered.biom -o 07_OTU_table/tufA/otus_to_remove.txt --to-tsv --table-type="OTU table"


# Remove only OTUs representing less than 1% of the reads:
# Requires all_otus.txt + otus_to_remove.txt and filter_otus_to_exclude.py
# You may want to increase the threshold for beta diversity analysis.
# get all otus.txt
biom convert -i 07_OTU_table/16S/OTU_table_2.biom -o 07_OTU_table/16S/OTU_table_2.txt --to-tsv --table-type="OTU table"
biom convert -i 07_OTU_table/18S/OTU_table_2.biom -o 07_OTU_table/18S/OTU_table_2.txt --to-tsv --table-type="OTU table"
biom convert -i 07_OTU_table/UPA/OTU_table_2.biom -o 07_OTU_table/UPA/OTU_table_2.txt --to-tsv --table-type="OTU table"
biom convert -i 07_OTU_table/tufA/OTU_table_2.biom -o 07_OTU_table/tufA/OTU_table_2.txt --to-tsv --table-type="OTU table"
# Filter real contamninats by 1% 
filter_otus_to_exclude.py -e 07_OTU_table/16S/otus_to_remove.txt -a 07_OTU_table/16S/OTU_table_2.txt -o 07_OTU_table/16S/otus_to_really_exclude.txt -t 1
filter_otus_to_exclude.py -e 07_OTU_table/18S/otus_to_remove.txt -a 07_OTU_table/18S/OTU_table_2.txt -o 07_OTU_table/18S/otus_to_really_exclude.txt -t 1
filter_otus_to_exclude.py -e 07_OTU_table/UPA/otus_to_remove.txt -a 07_OTU_table/UPA/OTU_table_2.txt -o 07_OTU_table/UPA/otus_to_really_exclude.txt -t 1
filter_otus_to_exclude.py -e 07_OTU_table/tufA/otus_to_remove.txt -a 07_OTU_table/tufA/OTU_table_2.txt -o 07_OTU_table/tufA/otus_to_really_exclude.txt -t 1
# Filter otus based on the new otus_to_really_remove table
filter_otus_from_otu_table.py -i 07_OTU_table/16S/OTU_table_2.biom -o 07_OTU_table/16S/OTU_table_3.biom -e 07_OTU_table/16S/otus_to_really_exclude.txt
filter_otus_from_otu_table.py -i 07_OTU_table/18S/OTU_table_2.biom -o 07_OTU_table/18S/OTU_table_3.biom -e 07_OTU_table/18S/otus_to_really_exclude.txt
filter_otus_from_otu_table.py -i 07_OTU_table/UPA/OTU_table_2.biom -o 07_OTU_table/UPA/OTU_table_3.biom -e 07_OTU_table/UPA/otus_to_really_exclude.txt
filter_otus_from_otu_table.py -i 07_OTU_table/tufA/OTU_table_2.biom -o 07_OTU_table/tufA/OTU_table_3.biom -e 07_OTU_table/tufA/otus_to_really_exclude.txt

# Print second statistics
biom summarize_table -i 07_OTU_table/16S/OTU_table_3.biom -o 07_OTU_table/16S/OTUs_stats2.txt
biom summarize_table -i 07_OTU_table/18S/OTU_table_3.biom -o 07_OTU_table/18S/OTUs_stats2.txt
biom summarize_table -i 07_OTU_table/UPA/OTU_table_3.biom -o 07_OTU_table/UPA/OTUs_stats2.txt
biom summarize_table -i 07_OTU_table/tufA/OTU_table_3.biom -o 07_OTU_table/tufA/OTUs_stats2.txt

echo "Done!"



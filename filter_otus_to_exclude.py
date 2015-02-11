#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Delete OTUs from the contaminant OTUs list that are very frequent.
This is done to remove only the contaminants that are found in the water and in the reagents,
I still want to keep boring algae OTUs that ended up in the control - this is a problem too, 
but takig them out of the analysis means to exclude the most common boring algae.
USAGE: filter_otus_to_exclude.py -e otus_to_exclude.txt -a all_otus.txt -t 2

@author: VanessaRM
09 - Feb - 2014

"""
from argparse import ArgumentParser
from itertools import izip

parser = ArgumentParser()
parser.add_argument('-e', '--otus_to_exclude', help='The path to the input otus_to_exclude.txt', required=True)
parser.add_argument('-a', '--all_otus', help='The path to the all_otus.txt', required=True)
parser.add_argument('-o', '--output', help='The path and name of the output file', required=True)
parser.add_argument('-t', '--threshold_perc', help='The threshold percentage to keep or exclude the contamiant OTU', required=True)

args = parser.parse_args()
input_otus_to_exclude = args.otus_to_exclude
input_all_otus = args.all_otus
output_fp = args.output
threshold_perc = int(args.threshold_perc)

# Store all otus counts:
with open(input_all_otus) as file:
  all_otus = file.readlines()[2:]

dict_all_otus = {}
total_number_of_reads = []
for i in all_otus:
    i_split = i.split()
    i_ID = i_split[0]
    i_counts = i_split[1:]
    count_reads_per_otu = []
    for f in i_counts:
        f = int(float(f))
        count_reads_per_otu.append(f)
        reads_per_otu = sum(count_reads_per_otu)
    total_number_of_reads.append (reads_per_otu)
    dict_all_otus['%s' %i_ID] = [reads_per_otu]

#calculate the threshold to exclude OTUs in number of reads:
total_number_of_reads = sum(total_number_of_reads)
threshold_reads = total_number_of_reads * threshold_perc / 100

# OTUs from controls:
otus_IDs_to_exclude = []
def contaminat_otus_finder (otus_to_exclude):
    "get the otus ID of the otus_to_exclude file"
    with open(otus_to_exclude) as f:
        next(f) #skip the two first lines (headers)
        next(f)
        for otu in f:
            otu_split = otu.split('\t')
            otu_id = str(otu_split[0])
            otus_IDs_to_exclude.append(otu_id)
    return otus_IDs_to_exclude


# function to iterate over all otus, find the matching otus and calculate whether it should be kept or not.
def finder (all_ot, single_contam_otu):
    #print single_contam_otu
    counts = all_ot[single_contam_otu]
    counts = int(counts[0])
    if counts < threshold_reads:
        return "EXCLUDE"

    else:
        return "KEEP"


otus_to_really_exclude = []
def main (all_ot, contaminant_ot):
    for b in contaminant_ot:
        c = finder (all_ot, b)
        if c == "EXCLUDE":
            otus_to_really_exclude.append(b)
          
contaminat_otus_finder (input_otus_to_exclude)

if __name__ == '__main__':
    main(dict_all_otus, otus_IDs_to_exclude)

# Save file
count = 0
savefile = open(output_fp, "w")
for lines in otus_to_really_exclude:
    savefile.write("%s\n" %lines)
    count += 1
savefile.close()

print ""
print "Filtering done:"
print "Saved %i otus in the %s file." % (count, output_fp)


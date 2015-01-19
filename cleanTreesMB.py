# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 10:36:06 2014
Summarize BEAST trees by clade probabilities

"""
#import packages
import re
 
#inputs
file_path = "Combined_trees_Runs02_03_04_subsample_10000_final.tre"

trees_only = []
clean_trees = []
upper_text = []

dirty_trees = open(file_path, "r")

# Separate trees and upper text:
for line in dirty_trees:
    if line.find ("tree ") != -1:
        trees_only.append(line)
    else:
        upper_text.append (line)

#deal with the last "end"
del upper_text[-1] 
end = ["End;"]

for lines in trees_only:
    cleaning = re.sub ('\[.*?\]',"",lines) # removes everything inside brackets
    clean_trees.append (cleaning)

my_trees = upper_text + clean_trees + end

savefile = open("my_10000trees.tre", "w")
counter = 1
for lines in my_trees:
    print (counter)
    savefile.write("%s" % lines)
    counter += 1 

print ("done")


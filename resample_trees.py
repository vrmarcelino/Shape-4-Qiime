'''
Resample random trees from a nexus tree file
20-10-2013
'''

import re
from random import randint

clean_trees = []
trees_only = []
upper_text = []


#inputs
print ("")
file_path = raw_input ("Type the tree file path/name:")
print ("")
wanted_trees = raw_input ("How many trees do you wish?")
wanted_trees_num = int(wanted_trees)
print ("")


# Open trees file as plain text and delete the unwanted stuff
dirty_trees = open(file_path, "r")
for lines in dirty_trees:
    cleaning = re.sub ('\[.*?\]',"",lines) # removes everything inside brackets
    clean_trees.append (cleaning)


# Put all trees in a list
for line in clean_trees:
    if line.find ("STATE") != -1:
        trees_only.append(line)
    else:
        upper_text.append (line)

del upper_text[-1] # delete the last "End;"


# Delete random trees
remaining_trees = (len(trees_only))
while remaining_trees > wanted_trees_num:
    random_number = randint(0,((remaining_trees)-1)) 
    del trees_only[(random_number)]
    remaining_trees = (len(trees_only))


# Put lists together and save the file
end = ["End;"]
my_trees = upper_text + trees_only + end

savefile = open("my_trees.tre", "w")
for lines in my_trees:
    savefile.write("%s" % lines)


print ("Done. %i trees saved as 'my_trees.tre'" %(remaining_trees))


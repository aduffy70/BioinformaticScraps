#! /usr/bin/python

# comparelists.py
# Takes any number of files with lists of gene names and reports which
# genes are in each list and every combination of lists. Suitable for
# generating data to display on a venn diagram.

import sys

master_list = {} # A dictionary with key=gene name and value= a list showing whether the gene name is present in each file (1=True, 0=False)

list_count = len(sys.argv) - 1
for x in range(0,list_count):
    with open(sys.argv[x + 1]) as file_x:
        for line in file_x:
            if line[0] != "#":
                if line in master_list:
                    #update the dictionary to show that this gene is present in this file
                    master_list[line][x] = 1
                else:
                    #Add a blank entry of the correct length to the dictionary
                    master_list[line] = [0] * list_count
                    #Update the dictionary to show that this gene is present in this file
                    master_list[line][x] = 1
venn_data = {}
for gene in master_list.keys():
    if str(master_list[gene]) in venn_data:
        venn_data[str(master_list[gene])] += 1
    else:
        venn_data[str(master_list[gene])] = 1
print "Present in files", "Count_of_genes"
for combo in sorted(venn_data.keys()):
    print combo, venn_data[combo]
for combo in sorted(venn_data.keys()):
    print "\n\n\n\n", combo
    for gene in master_list.keys():
        if str(master_list[gene]) == combo:
            print gene,

#! /usr/bin/python

#compare2lists.py
# Takes 2 files with lists of gene names and reports which genes are in
# file1, which are in file2, and which are in both

import sys

list1 = {}

with open(sys.argv[1]) as file1:
    # make a dictionary of the contents of file 1
    for line in file1:
        if line[0] != "#":
            list1[line] = [1,0]
with open(sys.argv[2]) as file2:
    #check whether the item is in the dictionary. If not, add it. If so, update it.
    for line in file2:
        if line[0] != "#":
            if line in list1:
                list1[line][1] = 1
            else:
                list1[line] = [0,1]
in_both = 0
in_list1 = 0
in_list2 = 0
print "In_both_lists"
for gene in list1.keys():
    if list1[gene] == [1,1]:
        in_both += 1
        print gene,
print "\n\n\nIn_list_1"
for gene in list1.keys():
    if list1[gene] == [1,0]:
        in_list1 += 1
        print gene,
print "\n\n\nIn_list_2"
for gene in list1.keys():
    if list1[gene] == [0,1]:
        in_list2 += 1
        print gene,
print "\n\n"
print "Both, list1, list2"
print in_both, in_list1, in_list2

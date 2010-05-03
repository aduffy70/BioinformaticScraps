#! /usr/bin/env python
# fastafilter.py
# Filter a fasta file to include or exclude a list of sequences.
# Prints to STDOUT, so redirect to file.
# Usage: fastafilter.py <sequence or quality fasta file> <list of sequence identifiers> "yes" || "no"
# Aaron M Duffy aduffy70{at}gmail.com
# April 2010


# import modules
from sys import argv  # gives us a list of command line arguments
from Bio import SeqIO  # biopython tools for reading/parsing sequence files

# get file names
seqFileName = argv[1]
listFileName = argv[2]
# keep or exclude the list; default to keep
if (len(argv) > 3):
    keep = argv[3]
else:
    keep = "yes"

# get and store the list of sequence id's
idList = []
listFileHandle = open(listFileName, "r")
for line in listFileHandle:
    line = line.rstrip()  #remove the line endings
    idList.append(line)
listFileHandle.close()

# index the fasta file (acts like a dictionary but without the memory contraints)
# the sequence id is the dictionary key
seqRecords = SeqIO.index(seqFileName, "fasta")
if (keep == "yes"):
    for id in idList:
        print ">" + id
        print seqRecords[id].seq
else:
    for id in seqRecords.keys():
        if not(id in idList):
            print">" + id
            print seqRecords[id].seq
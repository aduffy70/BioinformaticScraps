#! /usr/bin/env python
# fastalengths.py
# Return the name and length of each sequence in a fastafile
# Prints to STDOUT, so redirect to file.
# Usage: fastalengths.py <sequence or quality fasta file>
# Aaron M Duffy aduffy70{at}gmail.com
# May 2010


# import modules
from sys import argv  # gives us a list of command line arguments
from Bio import SeqIO  # biopython tools for reading/parsing sequence files

# get file names
seqFileName = argv[1]

# index the fasta file (acts like a dictionary but without the memory contraints)
# the sequence id is the dictionary key
seqRecords = SeqIO.index(seqFileName, "fasta")
print "#name,length"
for seqKey in seqRecords.keys():
    length = len(seqRecords[seqKey].seq)
    name = seqRecords[seqKey].id
    print "%s,%s" % (name, length)


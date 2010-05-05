#! /usr/bin/env python
# fastafilter.py
# Extracts a single entry from a fastafile
# Usage: fastafilter.py <sequence or quality fasta file> <id of the entry to extract>
# Aaron M Duffy aduffy70{at}gmail.com
# May 2010


# import modules
from sys import argv  # gives us a list of command line arguments
from Bio import SeqIO  # biopython tools for reading/parsing sequence files

# get file names
seqFileName = argv[1]
extractEntry = argv[2]

# index the fasta file (acts like a dictionary but without the memory contraints)
# the sequence id is the dictionary key
seqRecords = SeqIO.index(seqFileName, "fasta")
if seqRecords.has_key(extractEntry):
    print ">" + extractEntry
    print seqRecords[extractEntry].seq

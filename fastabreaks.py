#! /usr/bin/env python
# fastabreaks.py
# Make sure that a fasta file has only 2 lines per record: 1 for the id and one for the sequence.  This avoids issues that extra line breaks can cause for some scripts.
# Prints to STDOUT, so redirect to file.
# Usage: fastabreaks.py <sequence or quality fasta file>
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
for seqKey in seqRecords.keys():
    print ">" + seqRecords[seqKey].id
    print seqRecords[seqKey].seq

#! /usr/bin/env python
# fastaglance.py
# Return some basic stats about a fastafile
# Prints to STDOUT, so redirect to file.
# Usage: fastaglance.py <sequence or quality fasta file>
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
recordCount = len(seqRecords)
print "# of records:", recordCount, " ",
longestSeq = ""
longestSeqLength = 0
shortestSeq = ""
shortestSeqLength = 100000000
sumSeqLength = 0
for seqKey in seqRecords.keys():
    length = len(seqRecords[seqKey].seq)
    sumSeqLength += length
    if (length <= shortestSeqLength):
        shortestSeqLength = length
        shortestSeq = seqRecords[seqKey].id
    if (length >= longestSeqLength):
        longestSeqLength = length
        longestSeq = seqRecords[seqKey].id
if (recordCount):
    print "Mean length:", sumSeqLength / recordCount
    print "Longest:", longestSeqLength, longestSeq
    print "Shortest:", shortestSeqLength, shortestSeq
else:
    print

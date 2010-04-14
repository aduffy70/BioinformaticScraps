#! /usr/bin/env python
# fastadiff.py
# Returns a list of sequences found in file1, but not file2 & vice versa
# Aaron M Duffy aduffy70{at}gmail.com
# April 2010


from Bio import SeqIO
from sys import argv
file1Name = argv[1]
file2Name = argv[2]

seq1Records = SeqIO.index(file1Name, "fasta")
seq2Records = SeqIO.index(file2Name, "fasta")
output = "In %(file2)s but not %(file1)s:" % {'file1': file1Name, 'file2': file2Name}
print output
for key in seq2Records.keys():
    if not(seq1Records.has_key(key)):
        print key
output = "In %(file1)s but not %(file2)s:" % {'file1': file1Name, 'file2': file2Name}
print "\n" + output
for key in seq1Records.keys():
    if not(seq2Records.has_key(key)):
        print key

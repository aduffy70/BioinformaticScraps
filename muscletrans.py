#! /usr/bin/env python
# muscletranny.py
# Given a fasta file of dna sequences, performs a translation alignment using muscle with default settings.
# Requires muscle in $PATH
# Uses the bacterial translation table
# Usage: muscletranny.py <file with multiple dna sequences in fasta format>
# Aaron M Duffy aduffy70{at}gmail.com
# May 2010


# import modules
from sys import argv  # gives us a list of command line arguments
from Bio import SeqIO  # biopython tools for reading/parsing sequence files
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
from Bio.Seq import SeqRecord
import subprocess
import sys

# Open the DNA sequence file
dnaSeqFile = open(argv[1], 'r')

# Read the fasta sequences into a list
dnaSeqDict = SeqIO.to_dict(SeqIO.parse(dnaSeqFile, "fasta"))

# Translate the sequences
aaSeqRecords = []
for key in dnaSeqDict:
    aaSeq = SeqRecord(dnaSeqDict[key].seq.translate(table=11), id=key)
    aaSeqRecords.append(aaSeq)
dnaSeqFile.close()
# temporarily write them to a file
#tempaafile = open('tempaa.txt', 'w')
#for aaSeq in aaSeqRecords:
#    print >>tempaafile, '>' + aaSeq.id + '\n' + aaSeq.seq
#tempaafile.close()
#dnaSeqFile.close()

# Align the aa sequences
commandLine = str(MuscleCommandline())
childProcess = subprocess.Popen(commandLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell=(sys.platform!="win32")) # note - it is important not to pipe stderr or it will hang
SeqIO.write(aaSeqRecords, childProcess.stdin, "fasta")
childProcess.stdin.close()
alignment = AlignIO.read(childProcess.stdout, "fasta")
print alignment[1].seq

# Convert the aa alignment back into dna
for taxon in alignment:
    pass


"""
# Align the sequences (dna now, need to change to aa)
commandLine = str(MuscleCommandline())
childProcess = subprocess.Popen(commandLine, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=(sys.platform!="win32"))
SeqIO.write(dnaSeqRecords, childProcess.stdin, "fasta")
childProcess.stdin.close()
alignment = AlignIO.read(childProcess.stdout, "fasta")
print alignment
"""


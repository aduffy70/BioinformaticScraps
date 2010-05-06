#! /usr/bin/env python
# msatblast.py
# For each blast hit in tabular blast output... Grabs two sequences from a blast hit, makes sure they are in the orientations they were in the blast output and aligns them using clustalw.
# Usage: msatblast.py <blast output file(tabular)> <taxon1 contig fasta file> <taxon2 contig fasta file>
# Note - I opened, manipulated, and saved the blast output file as a .csv so the taxon names have ""'s around them.  If this script was combined with conservedmsats.py, the blast output could be handled in xml format and make better use of Biopython tools. The 2 steps I performed manually (discarding blast hits < 100nt & removing duplicates could be done within the script and we wouldn't have to worry about the ""'s.
# Aaron M Duffy aduffy70{at}gmail.com
# May 2010

# import modules
import os  # tools for handling files and working on the command line
import subprocess  # tools for working on the command line
from sys import argv  # a list of command line arguments
import re  # tools for working with regular expressions
#from collections import defaultdict  # provides dictionary of lists
from Bio import SeqIO  # biopython tools for reading/parsing sequence files
from Bio.Align.Applications import ClustalwCommandline  # biopython tools for running clustalw


blastOutputFileName = argv[1]
contigFileName1 = argv[2]
contigFileName2 = argv[3]
blastOutputFileHandle = open(blastOutputFileName, 'r')
seqRecordsTaxon1 = SeqIO.index(contigFileName1, "fasta")
seqRecordsTaxon2 = SeqIO.index(contigFileName2, "fasta")
pattern = re.compile('"(\w+)"\s+"(\w+)"\s+\S+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s+(\d+)')
# group(1) is taxon1, group(2) is taxon2, group(3) is db start, group(4) is db stop
for line in blastOutputFileHandle:
    match = pattern.search(line)
    alignFastaFileName = 'z' + match.group(1) + '-' + match.group(2) + '.fasta'
    alignFastaFileHandle = open(alignFastaFileName, 'w')
    alignFastaFileHandle.write('>' + match.group(1) + '\n')
    alignFastaFileHandle.write(str(seqRecordsTaxon1[match.group(1)].seq) + '\n')
    alignFastaFileHandle.write('>' + match.group(2) + '\n')
    if (int(match.group(3)) < int(match.group(4))):
        alignFastaFileHandle.write(str(seqRecordsTaxon2[match.group(2)].seq) + '\n')
    else: # Need the reverse-complement
        alignFastaFileHandle.write(str(seqRecordsTaxon2[match.group(2)].seq.reverse_complement()) + '\n')
    alignFastaFileHandle.close()
    cmd = ClustalwCommandline(infile=alignFastaFileName)
    p = subprocess.Popen(str(cmd), shell=True)
    sts = os.waitpid(p.pid, 0)
    os.remove(alignFastaFileName)
    dndFileName = 'z' + match.group(1) + '-' + match.group(2) + '.dnd'
    os.remove(dndFileName)
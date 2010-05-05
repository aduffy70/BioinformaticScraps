#! /usr/bin/env python
# conservedmsats.py
# Identifies microsatellites potentially conserved between two taxa.
#    Parses output files from msatcommander and sorts the contigs based on the type
#    of repeat for two taxa.  For each repeat type found, blast the contigs from the
#    second taxon against a database of contigs from the first taxon.  Output is a
#    series of blast ouput files named by the repeat type.)
# Usage: conservedmsats.py <taxon1 msatcommander file> <taxon2 msatcommander file> <taxon1 contig fasta file> <taxon2 contig fasta file>
# Note - on OSX, msatcommander output files have mac type line endings which must be converted to unix type line endings before using this script
# Aaron M Duffy aduffy70{at}gmail.com
# May 2010

# import modules
import os  # tools for handling files and working on the command line
import subprocess  # tools for working on the command line
from sys import argv  # a list of command line arguments
import re  # tools for working with regular expressions
from collections import defaultdict  # provides dictionary of lists
from Bio import SeqIO  # biopython tools for reading/parsing sequence files
from Bio.Blast.Applications import NcbiblastnCommandline  # biopython tools for running local blastn

def SortContigsByType(mscFileName):
    # Set up a dictionary of lists (A list of contig names that have each repeat type)
    contigsByType = defaultdict(list)
    # get file names and open the files for reading
    mscFileHandle = open(mscFileName, 'r')

    # Define the regular expressions we want to match.
    pattern1 = re.compile('(\w+),.*,\((\w+)\)\^') # group(1) is the contig name, group(2) is the repeat type
    pattern2 = re.compile('complement of (\w+)') # group(1) is the reversed repeat type

    # We want to minimise the number of repeat types we are tracking, so for forward repeats we will use the repeat type (group(2) from pattern1) but for reversed repeats we will use the reversed repeat type (group(1) from pattern2).

    # Look for the regular expressions in each line of the file
    # Add the contig name to the list of contigs with that repeat type
    for line in mscFileHandle:
        match1 = pattern1.search(line)
        if match1: # The line is describing a repeat but we don't yet know if it is reversed
            match2 = pattern2.search(line)
            if match2: # It is a reversed repeat
                contigsByType[match2.group(1)].append(match1.group(1))
            else: # It is not a reversed repeat
                contigsByType[match1.group(2)].append(match1.group(1))
    mscFileHandle.close()
    return contigsByType

def FilterSpecificContigs(contigFileName, contigList, outputFileName):
    outputFileHandle = open(outputFileName, 'w')
    seqRecords = SeqIO.index(contigFileName, "fasta")
    for contig in contigList:
        outputFileHandle.write('>' + contig + '\n' + str(seqRecords[contig].seq) + '\n')
    outputFileHandle.close()



mscFileName1 = argv[1]
mscFileName2 = argv[2]
contigFileName1 = argv[3]
contigFileName2 = argv[4]
taxon1ContigsByType = SortContigsByType(mscFileName1)
taxon2ContigsByType = SortContigsByType(mscFileName2)
# For each repeat type, format a blast database using the sequences from one taxon and blast the sequences from the other taxon against it.
for type in taxon1ContigsByType:
    if taxon2ContigsByType.has_key(type):
        print type, len(taxon1ContigsByType[type]), len(taxon2ContigsByType[type])
        # Make the fasta sequence file for the blast database
        FilterSpecificContigs(contigFileName1, taxon1ContigsByType[type], 'tempblastdatabase')
        # Format the blast database
        cmd = 'formatdb -i tempblastdatabase -p F -n tempblastdatabase'
        p = subprocess.Popen(cmd, shell=True)
        sts = os.waitpid(p.pid, 0) #Makes this script wait til the database is ready
        # Make the fasta sequence file for the blast query
        FilterSpecificContigs(contigFileName2, taxon2ContigsByType[type], 'tempblastquery')
        # Do the blast search
        blastOutFileName = 'blastout' + type + '.txt'
        cmd = NcbiblastnCommandline(query='tempblastquery', db='tempblastdatabase', evalue=0.001, outfmt=6, out=blastOutFileName)
        p = subprocess.Popen(str(cmd), shell=True)
        sts = os.waitpid(p.pid, 0)
# Cleanup temporary files
os.remove('tempblastdatabase')
os.remove('tempblastquery')
os.remove('tempblastdatabase.nhr')
os.remove('tempblastdatabase.nin')
os.remove('tempblastdatabase.nsq')

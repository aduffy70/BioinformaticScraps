#! /usr/bin/env python
# sortbycsv.py
# Take a single file of fasta sequences and break them out into files based on a ccsv file.  Each row in the csv file contains a comma separated list of sequence names that should be in a file together.
# Example: For several cpGenomes we have a fasta file of all the gene sequences in that genome.  We want to have a file for each gene containing the sequence for that gene from each genome, but the names of the genes in these files are not consistent. We create a spreadsheet table with a column for each taxon and a row for each gene.  The cells contain the sequence names.  If we output that file as a csv, and concatenate all the taxon fasta files into a single fasta file, these can be used with this script to create the files by gene.
# Prints to STDOUT, so redirect to file.
# Usage: sortbycsv.py <fasta file of all genes for all taxa in any order> <csv file>
# Aaron M Duffy aduffy70{at}gmail.com
# June 2010

#import modules
from sys import argv  # gives us a list of command line arguments
from Bio import SeqIO  # biopython tools for reading/parsing sequence files


#Use the full fasta description line as the dictionary key (avoids problems with spaces)
def GetFullName(record):
    return record.description


fastaFile = open(argv[1],'r')
fastaRecords = SeqIO.to_dict(SeqIO.parse(fastaFile, 'fasta'), key_function=GetFullName)
csvFile = open(argv[2], 'r')

print "Read %s fasta records." % len(fastaRecords.keys())
print "CSV entries with no matching fasta record:"
for line in csvFile:
    fastaNames = line.split(',')
    if (fastaNames[0]):
        geneFile = open(fastaNames[0] + ".fasta",'w')
        for name in fastaNames[1:]:
            name = name.rstrip()
            if (name):
                if (fastaRecords.has_key(name.lstrip('>'))):
                    print >>geneFile, name
                    print >>geneFile, fastaRecords[name.lstrip('>')].seq
                else:
                    print name

        geneFile.close()
csvFile.close()
fastaFile.close()
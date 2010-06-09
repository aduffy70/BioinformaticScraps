#! /usr/bin/env python
# listgenes.py
#This script produces a non-redundant list of the names of the features in a genbank genome flat file
#Usage: listgenes.py <genbank flat file name> <optional additional genbank flat file names...>
#Aaron M. Duffy  aduffy70{at}gmail.com
#June 2010

from Bio import SeqIO # tools for parsing genbank files
from sys import argv  # a list of command line arguments
import re  # tools for working with regular expressions

for fileName in argv[1:]:
    #Read the genbank flat file
    gbFile = open(fileName, 'r')
    gbRecord = SeqIO.read(gbFile, 'genbank')
    gbFile.close()

    #Print the name of the genbank record
    print gbRecord.name

    #Print the name of each unique feature except the first one (it is summary data for the whole sequence)
    geneList = []
    for feature in gbRecord.features[1:]:
        #print "type:: ",
        #print feature.type
        if (not (feature.type == "misc_feature") and not(feature.type == "repeat_region")):
            if ("gene" in feature.qualifiers.keys()):
                #print "\tgene:: " + feature.qualifiers["gene"][0]
                if (feature.qualifiers["gene"][0] not in geneList):
                    geneList = geneList + [feature.qualifiers["gene"][0]]
                #else:
                    #print "duplicate"
            elif ("note" in feature.qualifiers.keys()):
                #print "\tnote:: " + feature.qualifiers["note"][0]
                if (feature.qualifiers["note"][0] not in geneList):
                    geneList = geneList + [feature.qualifiers["note"][0]]
                #else:
                    #print "duplicate"
            #else:
                #print "\tMissed it:: ",
                #print feature.qualifiers
    for gene in geneList:
        print gene
    print str(len(geneList)) + "\n"

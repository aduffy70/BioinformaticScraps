#! /usr/bin/env python
# listgenesPith.py
#This script produces a non-redundant fasta file of the features in a genbank genome flat file (Optimized for Pinus thunbergii and other genomes with the same annotation idiosyncrasies)
#Usage: listgenesPith.py <genbank flat file name> <taxon tag>
#Aaron M. Duffy  aduffy70{at}gmail.com
#June 2010

from Bio import SeqIO # tools for parsing genbank files
from Bio.Seq import Seq
from sys import argv  # a list of command line arguments
import re  # tools for working with regular expressions

def SliceSequence(fullSequence, location, strand):
    #forward strand
    if (strand == 1):
        return fullSequence[location.nofuzzy_start:location.nofuzzy_end]
    #reverse strand
    else:
        return fullSequence[location.nofuzzy_start:location.nofuzzy_end].reverse_complement()

#Read the genbank flat file
gbFile = open(argv[1], 'r')
gbRecord = SeqIO.read(gbFile, 'genbank')
gbFile.close()
taxonTag = argv[2]

#Print the name of the genbank record
print gbRecord.name

#Print the name of each CDS, tRNA, rRNA, or pseudogene feature
#Don't eliminate duplicates - the tRNA's are not annotated with anti-codons!
for feature in gbRecord.features[1:]:
    # Is it a pseudogene?
    if ((feature.type == "gene") and ("pseudo" in feature.qualifiers.keys())):
        print ">" + taxonTag + "_" + feature.qualifiers["gene"][0] + " (pseudogene)"
        print SliceSequence(gbRecord.seq, feature.location, feature.strand)
    # Is it a CDS, tRNA, or rRNA?
    elif ((feature.type == "CDS") or (feature.type == "tRNA") or (feature.type == "rRNA")):
        #Is there a /gene annotation?
        if ("gene" in feature.qualifiers.keys()):
            print ">" + taxonTag + "_" + feature.qualifiers["gene"][0]
        #Is there a /product annotation?
        elif ("product" in feature.qualifiers.keys()):
            print ">" + taxonTag + "_" + feature.qualifiers["product"][0]
        #Seriously?! Is there at least a /note annotation?
        elif ("note" in feature.qualifiers.keys()):
            print ">" + taxonTag + "_" + feature.qualifiers["note"][0]
        #I give up!
        else:
            print "******* ERROR *********"
        #We need to join multiple sequences
        if (len(feature.sub_features) > 0):
            exonSequences = []
            for subfeature in feature.sub_features:
                exonSequences = exonSequences + [SliceSequence(gbRecord.seq, subfeature.location, subfeature.strand)]
#                print "subfeature ", subfeature.location, " ", subfeature.strand
            joinedExons = Seq("")
            for exon in exonSequences:
                joinedExons = joinedExons + exon
            print joinedExons
        #Just one sequence fragment
        else:
            print SliceSequence(gbRecord.seq, feature.location, feature.strand)
#        print feature.location, " ", feature.strand
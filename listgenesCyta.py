#! /usr/bin/env python
# listgenesCyta.py
#This script produces a non-redundant fasta file of the features in a genbank genome flat file (Optimized for Cycas taitungensis and other genomes with the same annotation idiosyncrasies)
#Usage: listgenesCyta.py <genbank flat file name> <taxon tag>
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
#We can eliminate duplicates- the tRNA names include the anticodon
uniqueFeatures = []
for feature in gbRecord.features[1:]:
    duplicate = False
    # Is it a pseudogene?
    if ((feature.type == "gene") and ("pseudo" in feature.qualifiers.keys())):
        if not(feature.qualifiers["gene"][0] in uniqueFeatures):
            uniqueFeatures = uniqueFeatures + [feature.qualifiers["gene"][0]]
            print ">" + taxonTag + "_" + feature.qualifiers["gene"][0] + " (pseudogene)"
            print SliceSequence(gbRecord.seq, feature.location, feature.strand)
        else:
            duplicate = True
    # Is it a CDS, tRNA, or rRNA?
    elif ((feature.type == "CDS") or (feature.type == "tRNA") or (feature.type == "rRNA")):
        if (feature.type == "tRNA"):
            antiCodon = ""
            if ("codon_recognized" in feature.qualifiers.keys()):
                antiCodon = feature.qualifiers["codon_recognized"][0]
            if not((feature.qualifiers["gene"][0] + "-" + antiCodon) in uniqueFeatures):
                uniqueFeatures = uniqueFeatures + [feature.qualifiers["gene"][0] + "-" + antiCodon]
                print ">" + taxonTag + "_" + feature.qualifiers["gene"][0] + "-" + antiCodon
            else:
                duplicate = True
        #Is there a /gene annotation?
        elif ("gene" in feature.qualifiers.keys()):
            if not(feature.qualifiers["gene"][0] in uniqueFeatures):
                uniqueFeatures = uniqueFeatures + [feature.qualifiers["gene"][0]]
                print ">" + taxonTag + "_" + feature.qualifiers["gene"][0]
            else:
                duplicate = True
        #Is there a /note annotation?
        elif ("note" in feature.qualifiers.keys()):
            if not(feature.qualifiers["note"][0] in uniqueFeatures):
                uniqueFeatures = uniqueFeatures + [feature.qualifiers["note"][0]]
                print ">" + taxonTag + "_" + feature.qualifiers["note"][0]
            else:
                duplicate = True
        #Seriously?! Is there at least a /product annotation?
        elif ("product" in feature.qualifiers.keys()):
            if not(feature.qualifiers["product"][0] in uniqueFeatures):
                uniqueFeatures = uniqueFeatures + [feature.qualifiers["product"][0]]
                print ">" + taxonTag + "_" + feature.qualifiers["product"][0]
            else:
                duplicate = True
        #I give up!
        else:
            print "\n\n\n******* ERROR *********\n\n\n"
        if not duplicate:
            #We need to join multiple sequences
            if (len(feature.sub_features) > 0):
                exonSequences = []
                for subfeature in feature.sub_features:
                    exonSequences = exonSequences + [SliceSequence(gbRecord.seq, subfeature.location, subfeature.strand)]
#                    print "subfeature ", subfeature.location, " ", subfeature.strand
                joinedExons = Seq("")
                for exon in exonSequences:
                    joinedExons = joinedExons + exon
                print joinedExons
            #Just one sequence fragment
            else:
                print SliceSequence(gbRecord.seq, feature.location, feature.strand)
#            print feature.location, " ", feature.strand
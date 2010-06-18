#! /usr/bin/env python
# gb2tbl.py
#This script converts a genbank flat file to a features table suitable for use with Sequin.
#Usage gb2tbl.py <genbank flatfile name>
#Writes to standard output so redirect to a file if desired
#Aaron M. Duffy  aduffy70{at}gmail.com
#May 2010

from Bio import SeqIO # tools for parsing genbank files
from sys import argv  # a list of command line arguments
import re  # tools for working with regular expressions

#Read the genbank flat file
gbFile = open(argv[1], 'r')
gbRecord = SeqIO.read(gbFile, 'genbank')

#Print the header row
print ">Feature gb|%s|" % gbRecord.name

#Setup a pattern match to filter out "Geneious name:" lines
pattern = re.compile('Geneious name')

#Format and print each feature except the first one (it is summary data for the whole sequence)
for feature in gbRecord.features[1:]:
    if (len(feature.sub_features) > 0): # Handle features with no subfeatures
        firstSubFeature = True;
# It looks like this step of reversing the order of subfeatures for features on the reverse strand is not necessary.  When this code was included and we imported to Sequin, all the reverse strand features had their exons in the wrong order.  I'm commenting this out but leaving it in til it has actually been tested and confirmed to be unnecessary.
#        if (feature.strand == -1): # Reverse the order of subfeatures on the reverse strand
#            orderedSubfeatures = reversed(feature.sub_features)
#        else:
#            orderedSubfeatures = feature.sub_features
        orderedSubfeatures = feature.sub_features #
        for subfeature in orderedSubfeatures:
            if (subfeature.strand == -1): # reverse strand
                start = subfeature.location.nofuzzy_end
                stop = subfeature.location.nofuzzy_start + 1 # adjust for the python 0-index
            else: # forward strand
                start = subfeature.location.nofuzzy_start + 1 # adjust for the python 0-index
                stop = subfeature.location.nofuzzy_end
            if (firstSubFeature): # Only print the subfeature type for the first subfeature
                print "%s\t%s\t%s" % (start, stop, subfeature.type)
                firstSubFeature = False
            else:
                print "%s\t%s" % (start, stop)
            for key in subfeature.qualifiers.keys():
                if ((key != "codon_start") and (key != "transl_table") and (key != "translation")):
                    print "\t\t\t%s%t%s" % (key, subfeature.qualifiers[key][0])
        for key in feature.qualifiers.keys():
            if ((key != "codon_start") and (key != "transl_table") and (key != "translation") and (key != "db_xref")and (key != "modified_by") and (key != "created_by")):
                if not (pattern.search(feature.qualifiers[key][0])):
                    if (key == "protein_id"):
                        print "\t\t\t%s\tgb|%s|" % (key, feature.qualifiers[key][0])
                    else:
                        print "\t\t\t%s\t%s" % (key, feature.qualifiers[key][0])
    else: # handle features with subfeatures
        if (feature.strand == -1): # reverse strand
            start = feature.location.nofuzzy_end
            stop = feature.location.nofuzzy_start + 1 # adjust for the python 0-index
        else: # forward strand
            start = feature.location.nofuzzy_start + 1 # adjust for the python 0-index
            stop = feature.location.nofuzzy_end
        print "%s\t%s\t%s" % (start, stop, feature.type)
        for key in feature.qualifiers.keys():
            if ((key != "codon_start") and (key != "db_xref") and (key != "transl_table") and (key != "translation") and (key != "modified_by") and (key != "created_by")): #and not((key == "gene") and ((feature.type == "tRNA") or (feature.type == "rRNA") or (feature.type == "CDS")))):
                if not (pattern.search(feature.qualifiers[key][0])):
                    if (key == "protein_id"):
                        print "\t\t\t%s\tgb|%s|" % (key, feature.qualifiers[key][0])
                    else:
                        print "\t\t\t%s\t%s" % (key, feature.qualifiers[key][0])

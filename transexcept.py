#! /usr/bin/env python
#transexcept.py
#Quick and dirty script to transfer translation exception info from the miscellaneous feature to the CDS in a genbank flat file
#Usage: transexcept.py <gff-like file with just the translation exception lines> <feature table-like file with just the CDS sections>
#Writes to standard output so redirect to a file if desired
#Aaron M. Duffy  aduffy70{at}gmail.com
#June 2010

import re
from sys import argv  # a list of command line arguments
from collections import defaultdict # provides a dictionary of lists

gffFile = open(argv[1], 'r')
cdsFile = open(argv[2], 'r')


#Read the gff file data into a dictionary
gffPattern = re.compile('(\d+)\s+(\d+).+([\-\+]).+transl_except(\w)')
gffData = defaultdict(list)
for line in gffFile:
    match = gffPattern.search(line)
    if (match):
            gffData[match.group(1)] = [match.group(1), match.group(2), match.group(3), match.group(4)]

#Cycle through the cdsFile
numberLinePattern = re.compile('(\d+)\s+(\d+)')
otherLinePattern = re.compile('^\s+')
exceptions = []
for line in cdsFile:
    match = numberLinePattern.search(line)
    if (match):
        print line,
        if (int(match.group(1)) < int(match.group(2))):
            start = int(match.group(1))
            end = int(match.group(2))
        else:
            start = int(match.group(2))
            end = int(match.group(1))
        for key in gffData.keys():
            if ((int(key) >= start) and (int(key) <= end)):
                exceptions.append(gffData[key])
    else:
        match = otherLinePattern.search(line)
        if (match):
            if (len(exceptions) > 0):
                for exception in exceptions:
                    if (exception[2] == '+'):
                        print '\t\t\ttransl_except\t"(pos:' + exception[0] + '..' + exception[1] + ', aa:' + exception[3] + ')"'
                    else:
                        print '\t\t\ttransl_except\t"(pos:' + exception[1] + '..' + exception[0] + ', aa:' + exception[3] + ')"'
                exceptions = [] #clear the list of exceptions
            print line,

gffFile.close()
cdsFile.close()
#! /usr/bin/env python
# wilcoxparse.py
# Extracts just the p-values from the output of wilcox.py.  Output of this script is csv format: GeneName,dN_p-value,dS_p-value
# Use: wilcox.py <path to text file with wilcox.py output>
# Aaron M Duffy aduffy70{at}gmail.com
# July 2010

#import modules
from sys import argv  # gives us a list of command line arguments
import os # for handling file paths 
import glob # lets us handle each file in a directory of files
import re # regular expressions
import subprocess # for running shell commands

path = argv[1] # Path to text file holding wilcox.py output
inFile = open(path, 'r')
geneNames = []
dNpValues = []
dSpValues = []
dNdSType = ""
genePattern = re.compile('FilesbyGene/(\S+)\.csv')
dNPattern = re.compile('Fern_dN and Seedplant_dN')
dSPattern = re.compile('Fern_dS and Seedplant_dS')
pValuePattern = re.compile('p-value = (\S+)$')
for line in inFile:
    match = genePattern.search(line)
    if (match):
        geneNames.append(match.group(1))
    else:
        match = dNPattern.search(line)
        if (match):
            dNdSType = "dN"
        else:
            match = dSPattern.search(line)
            if (match):
                dNdSType = "dS"
            else:
                match = pValuePattern.search(line)
                if (match):
                    if (dNdSType == "dN"):
                        dNpValues.append(match.group(1))
                    elif (dNdSType == "dS"):
                        dSpValues.append(match.group(1))
                    else:
                        print "ERROR"
for x in range(len(geneNames)):
    print geneNames[x] + "," + dNpValues[x] + "," + dSpValues[x]


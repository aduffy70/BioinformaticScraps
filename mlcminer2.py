#! /usr/bin/env python
# mlcminer2.py
# Extracts dN and dS values from a directory of PAML output files (.mlc files) and writes them to a separate csv file for each gene.  It "mines" data from mlc files.  Unlike mlcminer.py, this script doesn't average the dN and dS values - it extracts the raw values so we can use them in statistical analyses.  It also ignores the outgroup values.
# Aaron M Duffy aduffy70{at}gmail.com
# July 2010

#import modules
from sys import argv  # gives us a list of command line arguments
import os # for handling file paths
import glob # lets us handle each file in a directory of files
import re # regular expressions
from collections import defaultdict # provides a dictionary of lists

path = argv[1] # Path to a folder holding multiple PAML output (.mlc) files


#Get the dN values
pattern = re.compile('^\s+\d+\.\.\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+$')
for fileName in glob.glob(os.path.join(path, '*.mlc')):
    geneName = fileName.split('/')[-1].split('.')[0]
    firstTime = True
    Previous = 0
    group = 0
    maxLength = 0
    groupNames = ["Outgroup", "Seedplants", "Monilophytes"]
    dNValues = defaultdict(list)
    dSValues = defaultdict(list)
    mlcFile = open(fileName, 'r')
    for line in mlcFile:
        match = pattern.search(line)
        if (match):
            dNdS = float(match.group(1))
            dN = float(match.group(2))
            dS = float(match.group(3))
            if (firstTime or (dNdS == Previous)):
                dNValues[groupNames[group]].append(dN)
                dSValues[groupNames[group]].append(dS)
                firstTime = False
                Previous = dNdS
            else:
                group+=1
                dNValues[groupNames[group]].append(dN)
                dSValues[groupNames[group]].append(dS)
                Previous = dNdS
    mlcFile.close()
    if (len(dNValues["Seedplants"]) > len(dNValues["Monilophytes"])):
       maxLength = len(dNValues["Seedplants"])
    else:
       maxLength = len(dNValues["Monilophytes"])
    outFile = open(geneName + '.csv', "w")
    print >>outFile, "Fern_dN,Seedplant_dN,Fern_dS,Seedplant_dS"
    for x in range(maxLength):
       try:
           SeeddN = dNValues["Seedplants"][x]
       except:
           SeeddN = ""
       try:
           FerndN = dNValues["Monilophytes"][x]
       except:
           FerndN = ""
       try:
           SeeddS = dSValues["Seedplants"][x]
       except:
           SeeddS = ""
       try:
           FerndS = dSValues["Monilophytes"][x]
       except:
           FerndS = ""
       print >>outFile, str(FerndN) + "," + str(SeeddN) + "," + str(FerndS) + "," + str(SeeddS)
    outFile.close()


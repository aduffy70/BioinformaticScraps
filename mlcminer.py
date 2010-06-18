#! /usr/bin/env python
# mlcminer.py
# Extracts dN/dS values from a directory of PAML output files (.mlc files).  It "mines" data from mlc files...
# 3 scripts in one.  Uncomment one section at a time to get dNdS values, dN values, and dS values
# Aaron M Duffy aduffy70{at}gmail.com
# June 2010

#import modules
from sys import argv  # gives us a list of command line arguments
import os # for handling file paths
import glob # lets us handle each file in a directory of files
import re # regular expressions
from collections import defaultdict # provides a dictionary of lists

path = argv[1]

"""
#Get the dNdS values
pattern = re.compile('w \(dN/dS\) for branches\:\s+(\S+)\s+(\S+)\s+(\S+)')
for fileName in glob.glob(os.path.join(path, '*.mlc')):
    geneName = fileName.split('/')[1].split('.')[0]
    mlcFile = open(fileName, 'r')
    for line in mlcFile:
        match = pattern.search(line)
        if (match):
            print '%s,%s,%s,%s' % (geneName, match.group(1), match.group(2), match.group(3))
    mlcFile.close()
"""

"""
#Get the dN values
pattern = re.compile('^\s+\d+\.\.\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+$')
for fileName in glob.glob(os.path.join(path, '*.mlc')):
    geneName = fileName.split('/')[1].split('.')[0]
    print geneName,
    firstTime = True
    Previous = 0
    dNSum = 0
    count = 0
    mlcFile = open(fileName, 'r')
    for line in mlcFile:
        match = pattern.search(line)
        if (match):
            dNdS = float(match.group(1))
            dN = float(match.group(2))
            dS = float(match.group(3))
            if (firstTime or (dNdS == Previous)):
                dNSum = dNSum + dN
                count+=1
                firstTime = False
                Previous = dNdS
            else:
                dNAverage = dNSum / count
                print dNAverage,
                count = 1
                dNSum = dN
                Previous = dNdS
    dNAverage = dNSum / count
    print dNAverage
"""

#Get the dS values
pattern = re.compile('^\s+\d+\.\.\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+$')
for fileName in glob.glob(os.path.join(path, '*.mlc')):
    geneName = fileName.split('/')[1].split('.')[0]
    print geneName,
    firstTime = True
    Previous = 0
    dSSum = 0
    count = 0
    mlcFile = open(fileName, 'r')
    for line in mlcFile:
        match = pattern.search(line)
        if (match):
            dNdS = float(match.group(1))
            dN = float(match.group(2))
            dS = float(match.group(3))
            if (firstTime or (dNdS == Previous)):
                dSSum = dSSum + dS
                count+=1
                firstTime = False
                Previous = dNdS
            else:
                dSAverage = dSSum / count
                print dSAverage,
                count = 1
                dSSum = dS
                Previous = dNdS
    dSAverage = dSSum / count
    print dSAverage
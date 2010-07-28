#! /usr/bin/env python
# wilcox.py
# A python wrapper for the wilcox.r script.  That script runs two wilcoxon rank tests (1 on the first 2 columns of a csv file and one on the 3rd & 4th columns).  This wrapper lets us run the R script on every csv file in a folder and save the results to one long text file.
# Use: wilcox.py <path to folder with csv files> <path to wilcox.r>
# Writes to stdout, so redirect to a file.  When redirected, the warnings won't get redirected, which makes parsing the output easier.
# Aaron M Duffy aduffy70{at}gmail.com
# July 2010

#import modules
from sys import argv  # gives us a list of command line arguments
import os # for handling file paths 
import glob # lets us handle each file in a directory of files
import re # regular expressions
import subprocess # for running shell commands

csvpath = argv[1] # Path to a folder holding multiple csv files with dN dS values (the output of the mlcminer2.py script)
rpath = argv[2] # Path to the wilcox.r script
templist = ""
for fileName in glob.glob(os.path.join(csvpath, '*.csv')):
    geneName = fileName.split('/')[-1].split('.')[0]
    print "Gene Name = " + geneName
    os.system('R -f ' + rpath + ' --vanilla --slave --silent ' + fileName)


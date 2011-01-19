#!/usr/bin/env python
"""
Return some basic stats about a blast xml file
Prints to STDOUT, so redirect to file.
Usage: blastglance.py <blast xml file>
Aaron M Duffy aduffy70{at}gmail.com
January 2011
"""

from sys import argv
from Bio.Blast import NCBIXML

def main(args):
    if (len(args) > 1):
        input_filename = args[1]
    else:
        print("Error: No input or output filename provided")
        return 1
    blast_records = read_blast_results_from_file(input_filename)
    no_hits_count = 0 #Tracks the queries with no hits
    query_count = len(blast_records)
    for blast_record in blast_records:
        if (len(blast_record.descriptions) == 0):
            no_hits_count += 1
    print ('Queries: %s  Queries with no hits: %s' % (query_count,
                                                      no_hits_count))


def read_blast_results_from_file(filename):
    """
    Reads xml blast output from a file and returns a list of blast records.
    """
    blast_filehandle = open(filename, 'r')
    blast_records = NCBIXML.parse(blast_filehandle)
    blast_records_list = list(blast_records)
    return blast_records_list

if __name__ == '__main__':
    main(argv)

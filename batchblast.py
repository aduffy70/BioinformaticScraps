#!/usr/bin/env python
"""
batchblast.py
Takes a fasta file as input, runs blastn on all of the sequences and saves
the results to an xml file.  Also reads the results back in as a list
of blast record objects so you can process, display, or summarize as needed.
Usage: batchblast.py <inputfile> <outputfilename>
Aaron M. Duffy aduffy70{at}gmail.com
Jan 2010
"""

from sys import argv
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def main(args):
    ##################################
    #Customize these values as needed:
    e_value_cutoff = 0.0001
    number_of_hits = 10
    ##################################

    if (len(args) > 2):
        input_filename = args[1]
        output_filename = args[2]
    else:
        print("Error: No input or output filename provided")
        return 1
    print("Status: Reading input file")
    fasta_string = fasta_file_to_string(input_filename)
    print("Status: Running blast job")
    blast_output = run_ncbi_blast(fasta_string, e_value_cutoff,
                                         number_of_hits)
    print("Status: Saving blast results to file")
    save_blast_results_to_file(blast_output, output_filename)
    print("Status: Reading blast results from file")
    blast_results = read_blast_results_from_file(output_filename)
    """
    Insert commands to display, summarize, or otherwise process the blast
    results here.
    """
    print("Status: Complete")
    return 0


def fasta_file_to_string(filename):
    """
    Returns the fasta file as a single fasta format string.
    """
    fasta_filehandle = open(filename, 'r')
    fasta_string = fasta_filehandle.read()
    fasta_filehandle.close()
    return fasta_string

def run_ncbi_blast(fasta_string, e_value_cutoff, number_of_hits):
    """
    Blasts fasta sequences against NCBI's nr database.  Returns a blast
    output handle.  We could try to catch exceptions for more graceful
    exit on errors, but we'd lose error messages from NCBI.
    """
    blast_output = NCBIWWW.qblast('blastn', 'nr', fasta_string,
                                      expect=e_value_cutoff,
                                      hitlist_size=number_of_hits)
    return blast_output

def save_blast_results_to_file(blast_output, filename):
    """
    Saves xml blast output to a file.
    """
    blast_filehandle = open(filename, 'w')
    blast_filehandle.write(blast_output.read())
    blast_filehandle.close()

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

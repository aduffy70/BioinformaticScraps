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
    """Main function"""
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

    ### START SECTION- Uncomment to run one big blast job ###
    #fasta_string = fasta_file_to_string(input_filename)
    #print("Status: Running blast job")
    #blast_output = run_ncbi_blast(fasta_string, e_value_cutoff,
    #                                     number_of_hits)
    #print("Status: Saving blast results to file")
    #save_blast_results_to_file([blast_output], output_filename)
    ### END SECTION ###

    ### START SECTION- Uncomment to run many sequential blast jobs ###
    seqrecords_list = fasta_file_to_seqrecord_list(input_filename)
    print("Status: Running blast job")
    blast_output_list = run_sequential_ncbi_blast(seqrecords_list,
                                                  e_value_cutoff,
                                                  number_of_hits)
    print("Status: Saving blast results to file")
    save_blast_results_to_file(blast_output_list, output_filename)
    ### END SECTION ###

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

def save_blast_results_to_file(blast_output_list, filename):
    """
    Saves list of xml blast outputs to a file.  This still works for a
    single blast output - just submit it as a list with one element.
    """
    blast_filehandle = open(filename, 'w')
    for blast_output in blast_output_list:
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

def fasta_file_to_seqrecord_list(filename):
    """
    Returns the fasta file as a list of Seq record objects.
    """
    seq_records_list = list(SeqIO.parse(filename, 'fasta'))
    return seq_records_list

def run_sequential_ncbi_blast(seqrecords_list, e_value_cutoff,
                              number_of_hits):
    """
    Converts sequence records to fasta format and blasts them one-at-a-time
    against NCBI's nr database.  Provides status updates ('.' for success,
    '!' for failure) as each blast search completes. Writes any failing
    sequences to the console in fasta format. Returns a list of blast
    output handles.
    """
    blast_output_list = []
    failure_fasta_string = ''
    blast_count = 0
    for seqrecord in seqrecords_list:
        fasta_string = seqrecord.format("fasta")
        blast_count += 1
        try:
            blast_output_list.append(run_ncbi_blast(fasta_string,
                                                    e_value_cutoff,
                                                    number_of_hits))
            print(blast_count)
        except:
            print(str(blast_count) + " FAIL")
            failure_fasta_string += fasta_string
    if (failure_fasta_string != ''):
        print('\nFailures:')
        print(failure_fasta_string)
    else:
        print('\nNo blast failures')
    return blast_output_list



if __name__ == '__main__':
    main(argv)

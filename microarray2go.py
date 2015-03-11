#!/usr/bin/env python
"""
microarray2go.py
Takes a csv file of microarray markers as input and returns a file showing
the top blastn hit of the marker 70mers against an EST database and the
top blastx hit of those ESTs against the nr protein database.  Those
protein sequences can then be used in a searches against GO databases.

Aaron M. Duffy
December 2010
"""

from sys import argv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import re

Entrez.email = "aduffy70@gmail.com"

def main(args):
    # Hardcoded values for the Trout Annotation project.
    e_value_cutoff1 = 0.0001
    e_value_cutoff2 = 0.000001
    number_of_hits1 = 1
    number_of_hits2 = 10
    taxon_list = ['Oncorhyncus mykiss', 'Salmo salar', 'Danio rerio']
    filter_terms = ['buffer', 'empty']
    nr_filter_terms = ['predicted', 'unnamed protein product',
                       'hypothetical protein']
    """
    START SECTION: Run EST blast once then comment out.
    """
    """
    # Get the 70mer descriptions and sequences from the csv file
    seq_records = csv2seq(args[1], 5, 6, 13, 5, filter_terms, True)
    print('Read 70mers from csv file')
    # Store the filtered, unique records as a fasta file
    SeqIO.write(seq_records.values(), 'output/markers.fasta',
                'fasta')
    print('Stored 70mers as fasta file')
    # Blast the 70mers against a subset of NCBI's EST database
    blast_output_handle = est_blast(seq_records, taxon_list,
                                         e_value_cutoff1,
                                         number_of_hits1)
    print('Blasted 70mers against EST database')
    # Save blast results to file and read that for further use.
    blast_records = store_and_parse_blast_results(
                                    'output/EST_blastout.xml',
                                    [blast_output_handle])
    print('Stored blast results to xml file')
    """
    """
    END SECTION
    """

    """
    START TEMP: read the EST blast results from file so we can test without reblasting.  This section is only needed when the section above is commented out.
    """

    blast_output_handle = open('output/EST_blastout.xml', 'r')
    blast_records = NCBIXML.parse(blast_output_handle)
    blast_records = list(blast_records)

    """
    END TEMP
    """

    """
    START SECTION: Run nr blast once then comment out.
    """

    summarize_blast_results(blast_records)
    # Get a list of the GI numbers in the blast output
    blast_hit_gi_list = get_blast_gi_list(blast_records)
    print('Generated list of EST GI numbers')
    # Remove duplicates
    blast_hit_gi_list = list(set(blast_hit_gi_list))
    print('Removed dupicate GI numbers')
    # Get the fasta sequences for the GI numbers from NCBI
    get_genbank_records(blast_hit_gi_list)
    print('Downloaded EST sequences from NCBI')
    est_records = SeqIO.parse('output/EST_sequences.fasta', 'fasta')
    print('Read EST sequences from file')
    # Blast the est_records against nr
    blast_output_handle_list = nr_blast(est_records, e_value_cutoff2,
                                   number_of_hits2)
    print('Blasted ESTs against nr database')
    # Save blast results to file and read that for further use.
    blast_records = store_and_parse_blast_results(
                                'output/nr_blastout.xml',
                                blast_output_handle_list)
    print('Stored blast results to xml file')
    summarize_blast_results(blast_records)

    """
    END SECTION
    """

    """
    START TEMP: read the blast results from file so we can test without reblasting.  This section is only needed when the section above is commented out.
    """
    """
    blast_output_handle = open('output/nr_blastout.xml', 'r')
    blast_records = NCBIXML.parse(blast_output_handle)
    blast_records = list(blast_records)
    """
    """
    END TEMP
    """

    # Get the best hit from each blast result
    best_hits = get_best_blast_result(blast_records, nr_filter_terms)
    print('nr blast queries: %s ' % len(blast_records))
    print('nr blast queries with a best hit: %s' % len(best_hits))
    # Generate the results table
    create_results_table('output/EST_blastout.xml', best_hits)

def create_results_table(est_blast_file, best_hits):
    """
    Generates an output file showing each 70mer, its top EST blastn hit, 
    and the ESTs best nr blastx hit.
    """
    # Dictionary to store the EST hit for each 70mer.
    seventymer_hits = {}
    # Parse the est blast results
    est_blast_output_handle = open(est_blast_file, 'r')
    est_blast_records = NCBIXML.parse(est_blast_output_handle)
    est_blast_records = list(est_blast_records)
    # Store the EST hit for each 70mer or NA if none
    for blast_record in est_blast_records:
        if (len(blast_record.descriptions) > 0):
            for alignment in blast_record.alignments:
                seventymer_hits[blast_record.query] = alignment.hit_id
        else:
            seventymer_hits[blast_record.query] = 'NA'
    output_file = open('output/combined_blast_summary.csv', 'w')
    for key in seventymer_hits.keys():
        best_hit = 'NA'
        if (seventymer_hits[key] != 'NA'):
            seventymer_gi = seventymer_hits[key].split('|')[1]
            if seventymer_gi in best_hits.keys():
                best_hit = best_hits[seventymer_gi][1]
        output_file.write('%s,%s,%s\n' % (key, seventymer_hits[key], best_hit))
    output_file.close()

def get_best_blast_result(blast_records, filter_terms):
    """
    Returns the best blast hit for each query as a dictionary of query
    gi numbers (key) to tuples of the query accession description and the
    best hit accessions.  The 'best' hit is the first one that does not
    contain any term from the list of filter terms.
    """
    # Generate a regular expression string from the filter terms
    filter_string = '|'.join(filter_terms)
    filter_pattern = re.compile(filter_string)
    best_hits = {}
    found_best = False
    no_best_found = 0
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            if not(re.search(filter_pattern, alignment.hit_def.lower())):
                best_hits[blast_record.query.split('|')[1]] = (
                                blast_record.query, alignment.accession)
                found_best = True
                break
        if (not found_best):
            no_best_found += 1
    print('nr blast results with ALL hits filtered out: %s' % no_best_found)
    return best_hits


def store_and_parse_blast_results(file_name, blast_output_handle_list):
    """
    Stores blast results to file and then parses the file to a list we can
    reuse.  There may be a better way to do this, but we can only iterate
    the blast handle once, and saving to disk lets us avoid rerunning the
    same jobs on NCBI's servers while developing and debugging this script.
    """
    blast_output_file = open(file_name, 'w')
    for blast_output_handle in blast_output_handle_list:
        blast_output_file.write(blast_output_handle.read())
        blast_output_handle.close()
    blast_output_file.close()
    # Parse the blast results file and convert it to a list we can reuse
    blast_output_handle = open(file_name, 'r')
    blast_records = NCBIXML.parse(blast_output_handle)
    blast_records = list(blast_records)
    return blast_records

def get_genbank_records(accession_list):
    """
    Downloads fasta sequences from ncbi for a list of accessions (gi numbers)
    and writes them to file.
    """
    epost_results = Entrez.read(Entrez.epost("nucleotide",
                                             id=",".join(accession_list)))
    print("Posted")
    webenv = epost_results["WebEnv"]
    query_key = epost_results["QueryKey"]
    batch_size = 20
    out_file = open('output/EST_sequences.fasta', 'w')
    fetch_handle = Entrez.efetch(db='nucleotide', rettype='fasta',
                                 webenv=webenv, query_key=query_key)
    print("Fetched")
    data = fetch_handle.read()
    fetch_handle.close()
    out_file.write(data)
    out_file.close()

def get_blast_gi_list(blast_records):
    """
    Returns a list of the genbank GI numbers of the hits in a list of blast
    records.
    """
    blast_hit_gi_list = []
    for blast_record in blast_records:
        if (len(blast_record.descriptions) > 0):
            for alignment in blast_record.alignments:
                gi = alignment.hit_id.split('|')[1]
                blast_hit_gi_list.append(gi)
                print(gi)
    return blast_hit_gi_list

def csv2seq(csv_file, id_column=0, description_column = 1, seq_column=2,
            filter_column=0, filter_terms=[], filter_duplicates=True):
    """
    Returns a dictionary of Bio.SequenceRecords extracted from a csv file
    (tab delimited) of microarray marker descriptions.  Removes markers with
    duplicate names or with data in a filter column that matches any of a
    list of filter strings.
    """
    seq_records = {}
    file_handle = open(csv_file, 'r')
    for line in file_handle:
        csv_elements = line.rstrip().split('\t')
        if (len(csv_elements) < 13 ):
            print 'Not enough columns: %s' % csv_elements
        else:
            # Check against the list of filter terms
            if (csv_elements[filter_column].strip().lower()
                not in filter_terms):
                id = csv_elements[id_column].strip()
                # Check that it is not a duplicate
                if ((not filter_duplicates) or (id not in seq_records)):
                    seq = Seq(csv_elements[seq_column])
                    seq_records[id] = SeqRecord(seq, id=id,
                                                description=csv_elements[
                                                description_column])
    file_handle.close()
    return seq_records

def nr_blast(seq_records, expect_value=10, number_of_hits=50):
    """
    Uses a list of Bio.SeqRecords as the sequence queries and blasts
    NCBI's nr database.  Returns a list of xml format blast result handles.
    """
    # Create a fasta format string from a SeqRecord and blast it
    blast_results = []
    #Try only blasting against Danio rerio sequences
    for record in seq_records:
        fasta_string = record.format("fasta")
        try:
            print(fasta_string)
            blast_results.append(NCBIWWW.qblast("blastx", "nr", fasta_string,
                                   expect=expect_value,
                                   hitlist_size=number_of_hits))
        except:
            print('Failed: %s' % fasta_string)
    return blast_results

def est_blast(seq_records, organisms, expect_value=10, number_of_hits=50):
    """
    Uses a dictionary of Bio.SeqRecords as the sequence queries and blasts
    NCBI's EST database for a list of organisms.  Returns an xml format
    blast result handle.
    """
    # Create an entrez query string using the list of organisms
    entrez_query_string = '(%s[Organism]' % organisms[0]
    for taxon in organisms[1:]:
        entrez_query_string += ' OR %s[Organism]' % taxon
    entrez_query_string += ')'
    # Create a fasta format string from a list of SeqRecords and blast it
    fasta_string = ''
    for record in seq_records.values():
        fasta_string += record.format("fasta")
    blast_results = NCBIWWW.qblast("blastn", "est", fasta_string,
                         entrez_query=entrez_query_string,
                         expect=expect_value,
                         hitlist_size=number_of_hits)
    return blast_results

def summarize_blast_results(blast_records):
    no_hit_queries = []
    for blast_record in blast_records:
        if (len(blast_record.descriptions) == 0):
            no_hit_queries.append(blast_record.query)
    print('Total queries: %s\nQueries with no hits: %s' % (len(blast_records),
          len(no_hit_queries)))


if __name__ == '__main__':
    main(argv)

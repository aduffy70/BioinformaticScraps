#! /usr/bin/env python

# filter_structure_file_by_informative_loci.py
# Given a structure file (.str or .stru) remove all loci that are are not informative. To be informative a locus needs to be:
#  1. variable - there must be at least two genotypes in the population.
#  2. not autapomorphic - each genotype needs to be present in more than one sample.

# usage:
# filter_structure_file_by_informative_loci.py structure_file.stru output_stats_file_name > new_filtered_structure_file.stru

import os
import sys

filename = sys.argv[1]

locus_counts = {} # key = sample, value = [starting_loci_count, ending_loci_count]
total_starting_loci = 0 #Combined loci in all samples before filtering
total_remaining_loci = 0 # Combined loci in all samples after filtering
remaining_samples = 0 # Sample count after filtering
output_filename = sys.argv[2]
genotype_counts = [] # List with a dictionary for each locus. key=genotype, value=count of samples with that genotype
is_first_sample = True

#Get starting genotype counts for each locus
with open(filename, "r") as str_file:
    is_first_allele = True # 2 lines for each sample. We need both to get the genotype
    for line in str_file:
        elements = line.split()
        if is_first_allele:
            sample = elements[0]
            first_alleles = elements[1:]
            is_first_allele = False
        else:
            second_alleles = elements[1:]
            is_first_allele = True
            locus_count = 0
            for x in range(len(first_alleles)):
                genotype = str(min(first_alleles[x], second_alleles[x])) + str(max(first_alleles[x], second_alleles[x]))
                if is_first_sample:
                    genotype_counts.append({}) # Make a dictionary to hold the genotype counts for this locus
                if genotype != "-9-9":
                    if genotype in genotype_counts[x].keys():
                        genotype_counts[x][genotype] += 1
                    else:
                        genotype_counts[x][genotype] = 1
                    locus_count += 1
            is_first_sample = False
            locus_counts[sample] = [locus_count, 0]
    total_starting_loci = len(elements) - 1 #All lines have the same length so only count on the last one

#Now we know how many are present for each locus. Lets go back through the file again and only output loci that pass the filter

with open(filename, "r") as str_file:
    for line in str_file:
        elements = line.split()
        sample = elements[0]
        remaining_loci_count = 0 #Keep track of how many loci are present. If it stays zero, we can remove this sample entirely.
        updated_line = sample
        alleles = elements[1:]
        bad_over_3_genotypes = 0 # Count of loci with too many genotypes (why did ipyrad not filter these?)
        bad_1_genotype = 0 # Count of loci with no variation
        bad_2_genotypes = 0 # Count of loci with variation in only one sample
        good_2_genotypes = 0 # Count of loci with two genotypes each represented in >1 sample
        good_3_genotypes = 0 # Count of loci with three genotypes
        for x in range(len(alleles)):
            keep_locus = True
            if len(genotype_counts[x]) < 2: # No genotypes or 1 genotype. Omit locus.
                keep_locus = False
                bad_1_genotype += 1
                #print x, "BAD 0 or 1 genotypes", genotype_counts[x]
            elif len(genotype_counts[x]) == 3: #Three genotypes. Keep locus.
                keep_locus = True
                good_3_genotypes += 1
                #print x, "GOOD 3 genotypes", genotype_counts[x]
            elif len(genotype_counts[x]) == 2: # Two genotypes. Keep if both are present in at least two samples.
                for genotype in genotype_counts[x].keys():
                    if genotype_counts[x][genotype] < 2:
                        keep_locus = False
                if keep_locus:
                    good_2_genotypes += 1
                    #print x, "GOOD 2 genotypes", genotype_counts[x]
                else:
                    bad_2_genotypes += 1
                    #print x, "BAD 2 genotypes but <2 samples", genotype_counts[x]
            else:
                bad_over_3_genotypes += 1
                #print x, "BAD too many genotypes", genotype_counts[x]
            if keep_locus:
                updated_line += "\t" + alleles[x]
                if alleles[x] != "-9":
                    remaining_loci_count += 1
        locus_counts[sample][1] = remaining_loci_count
        if remaining_loci_count > 0:
            print updated_line
            remaining_samples += 1
    total_remaining_loci = len(updated_line.split()) - 1
    remaining_samples = remaining_samples / 2 #Because we counted each sample twice (2 lines per sample)

with open(output_filename, "w") as output_file:
    output_file.write(" ".join(sys.argv) + "\n\n")
    output_file.write("Sample,Loci_before,Loci_after\n")
    for sample in sorted(locus_counts.keys()):
        before = locus_counts[sample][0]
        after = locus_counts[sample][1]
        output_file.write(sample + "," + str(before) + "," + str(after) + "\n")
    output_file.write("\nTotal_loci_before," + str(total_starting_loci) + "\n")
    output_file.write("Total_loci_after," + str(total_remaining_loci) + "\n")
    output_file.write("Total_samples_after," + str(remaining_samples) + "\n")
    output_file.write("Bad-too_many_genotypes," + str(bad_over_3_genotypes) + "\n")
    output_file.write("Bad-no_variation," + str(bad_1_genotype) + "\n")
    output_file.write("Bad-autapomorphy," + str(bad_2_genotypes) + "\n")
    output_file.write("Good-2_genotypes," + str(good_2_genotypes) + "\n")
    output_file.write("Good-3_genotypes," + str(good_3_genotypes) + "\n")

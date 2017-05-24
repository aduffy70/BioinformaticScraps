#! /usr/bin/python

# extract_hplc_data.py
# Aaron M. Duffy
# 2017-02-16

# Extracts Peak number, R.Time, and Area values from hplc output files in the current folder in
# a specified order (by tube).
# Expects hplc output files to be named ending with .txt but gets the names from the Sample Name field in the file--not from the file name.
# If you have more than one sample with the same sample name (Two s1 standards, for example) you need to manually edit the files to make the names unique (s1, s1a, etc.).

# usage extract_hplc_data.py order_string > output_file.csv
# The order string consists of comma separated sample names in the order you want them returned (NO SPACES and NO DUPLICATE NAMES)
# To match our "normal" order:
# 2,10,18,3,11,19,4,12,20,5,13,21,1,9,17,6,14,22,7,15,s1,s1a,s2,s2a,s3,s3a,s4,s4a,s5,s5a

import sys
import glob # For looking at filenames

order_string = sys.argv[1]

#print(order_string)
hplc_dict = {} # key=Sample name, value=list of [Peak,R.Time,Area]

# Read data from all properly named files in the folder into a dictionary
for filename in glob.iglob('*.txt'):
    with open(filename) as hplc_file:
        is_peak_table = False
        for line in hplc_file:
            if "Sample Name" in line:
                name_elements = line.split()
                if len(name_elements) == 3:
                    sample_name = name_elements[-1]
                    if sample_name not in hplc_dict:
                        hplc_dict[sample_name] = []
                    else:
                        sys.stderr.write("Error: File has a duplicate sample name: " + filename + "\n")
                else:
                    sys.stderr.write("Warning: File doesn't include a valid sample name: " + filename + "\n")
            elif "Peak#" in line:
                is_peak_table = True
            else:
                if is_peak_table:
                    if len(line.split()) < 2:
                        is_peak_table = False
                    else:
                        line_elements = line.split()
                        hplc_dict[sample_name].append([line_elements[0], line_elements[1], line_elements[4], filename])

# Get the order of the files and print info in that order
sample_order = order_string.split(',')
order_count = 1
print("Order,Sample_name,Peak,R.Time,Area,File")
for sample in sample_order:
    if sample in hplc_dict:
        for peak in hplc_dict[sample]:
            peak_number = peak[0]
            r_time = peak[1]
            area = peak[2]
            file_string = peak[3]
            print("%s,%s,%s,%s,%s,%s" % (order_count, sample, peak_number, r_time, area, file_string))
    else:
        sys.stderr.write("Error: Sample name was not found: " + sample + "\n")
    order_count += 1

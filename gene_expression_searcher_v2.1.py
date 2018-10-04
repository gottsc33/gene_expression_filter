#!/usr/bin/env python

#Version 2.1 - additional filtering for comparisons

#This python script reads in a list of transcript sequence identifiers that are
#of interest to you (a 1 column text file with a unique ID number one each line
#example AT1G05489 or MD17G099443). Usefully for shrinking a dataset to a size
#that is easier to manage in identifying any significant differences in
#expression.

#Version 2.1 upgraded a single line of code that could potential create slowed processing speeds due to a large input file being left open


#Files for reading including gene_exp.diff, isoform_exp.diff or cds_exp.diff.

#Note all search queries terms must be in same format as
#displayed in the cuffdiff output files.

#Code by Chris Gottschalk - 8/24/17 - email questions: gottsc33@msu.edu

import sys
import os

#loading in target gene file

in_name = raw_input("\n Hello user, please type in the $PATH to your target transcript sequence ID file and press 'Enter':")
in_gene = open(in_name, "r")
queries = []
with in_gene as query_list:
    for query in query_list:
        queries.append(query.strip('\n'))
in_gene.close()
#print queries #optional print list of query transcript sequence ID numbers

#reading in cuffdiff output

matches = []
cuffdiff_in = raw_input("\n Please type in the $PATH to your cuffdiff output file and press 'Enter':")
exp_data = open(cuffdiff_in, "r")

for line in exp_data:
    current_line = line.split()
    for query in queries:
        if query in current_line:
            matches.append(line)
            #print matches
            break

#writing an output file that only contains information as origionally written by
#cuffdiff it will only contain lines with the unique transcript sequence ID
output = open('gene_expression_searcher.temp', 'w+')
for line in matches:
    output.write(line)
output.close()
exp_data.close()

#Setting additional filtering for select pairwise comparisons
comparison_matches = []
#User must define the comparison list below are the two need for the Gala & HC Transcriptomes
#comparison_list = ['NTFB\tTFB','NT15\tT15','NT35\tT35','NT50\tT50','NT70\tT70'] #honeycrisp
comparison_list = ['GA1DAT\tC1DAT','GA2DAT\tC2DAT','GA5DAT\tC5DAT','GA15DAT\tC15DAT'] #gala
temp_in = open('gene_expression_searcher.temp', 'r')

with temp_in as non_compare:
    hits = non_compare.readlines()

for line in hits:
    for comparison in comparison_list:
        if comparison in line:
            comparison_matches.append(line)
            break

temp_in.close()

output2 = open('gene_expression_filter.diff', 'w+')
for line in comparison_matches:
    output2.write(line)
output2.close()

#remove the temp file
os.remove('gene_expression_searcher.temp')

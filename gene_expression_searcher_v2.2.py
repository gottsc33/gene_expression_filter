#!/usr/bin/env python

#Version 2.2 - additional filtering for comparisons

#This python script reads in a list of transcript sequence identifiers that are
#of interest to you (a 1 column text file with a unique ID number one each line
#example AT1G05489 or MD17G099443). Usefully for shrinking a dataset to a size
#that is easier to manage in identifying any significant differences in
#expression.

#Version 2.2 minor upgrads to reports and file inputing 8/14/18
#Version 2.1 upgraded a single line of code that could potential create slowed processing speeds due to a large input file being left open 4/26/18


#Files for reading including gene_exp.diff, isoform_exp.diff or cds_exp.diff.

#Note all search queries terms must be in same format as
#displayed in the cuffdiff output files.

#Code by Chris Gottschalk - 8/24/17 - email questions: gottsc33@msu.edu

import sys
import os

#loading in target gene file

in_name = raw_input("\n Hello user, please type in the $PATH to your target transcript sequence ID file and press 'Enter': ")
assert os.path.exists(in_name), "I did not find the file at, "+str(in_name)
in_gene = open(in_name, "r")
print("Your input file has been found and is loading...")

queries = []
with in_gene as query_list:
    for query in query_list:
        queries.append(query.strip('\n'))
in_gene.close()
#print queries #optional print list of query transcript sequence ID numbers

print("Input sequence IDs indexed...")

#reading in cuffdiff output

matches = []
cuffdiff_in = raw_input("\n Please type in the $PATH to your cuffdiff output file and press 'Enter':")
assert os.path.exists(cuffdiff_in), "I did not find the file at, "+str(cuffdiff_in)
exp_data = open(cuffdiff_in, "r")
print("Your cuffdiff file has been found and is loading...")

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
print("writing temp file of matches...")
for line in matches:
    output.write(line)
output.close()
exp_data.close()

#Setting additional filtering for select pairwise comparisons
comparison_matches = []
#User must define the comparison list below are the two need for the Gala & HC Transcriptomes
#comparison_list = ['NTFB\tTFB','NT15\tT15','NT35\tT35','NT50\tT50','NT70\tT70'] #honeycrisp
#comparison_list = ['GA1DAT\tC1DAT','GA2DAT\tC2DAT','GA5DAT\tC5DAT','GA15DAT\tC15DAT'] #gala
comparison_list = raw_input("\n Please input all comparisons your would like to retrieve as written in cuffdiff ouput. example input format (include '' and regex tabs) - 'Control\tTreatment1, 'Control\tTreatment2' - press 'Enter when complete': ")
print("Thank you for entering comparisons, reading in temp file now...")
temp_in = open('gene_expression_searcher.temp', 'r')

with temp_in as non_compare:
    hits = non_compare.readlines()

print("Identifying target transcripts/genes and coorisponding comparisons...")
for line in hits:
    for comparison in comparison_list:
        if comparison in line:
            comparison_matches.append(line)
            break

temp_in.close()

print("Writing output file and removing temp files...")
output2 = open('gene_expression_filter.diff', 'w+')
for line in comparison_matches:
    output2.write(line)
output2.close()

#remove the temp file
os.remove('gene_expression_searcher.temp')

print("Thank you for choosing Gene Expression Searcher for all your cuffdiff ouput needs! Happy transcriptome analyzing!")

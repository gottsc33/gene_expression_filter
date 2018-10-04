#!/usr/bin/env python

#This python script reads in a list of expression data and extracts only
#significantly different expressed genes, loci or cds.

#Code by Chris Gottschalk - 8/25/17 - email questions: gottsc33@msu.edu

import sys

#reading in cuffdiff output

matches = []
significant = ['yes',]
cuffdiff_in = raw_input("\n Please type in the $PATH to your filtered cuffdiff output file generated from gene_expression_searcher.py and press 'Enter':")
exp_data = open(cuffdiff_in, "r")
with exp_data as f:
    f = f.readlines()

for line in f:
    for phrase in significant:
        if phrase in line:
            matches.append(line)
            break
exp_data.close()

#writing an output file that only contains information as origionally written by
#cuffdiff it will only contain lines with the unique transcript sequence ID
output = open('significant_searcher.out', 'w+')
for line in matches:
    output.write(line)
output.close()

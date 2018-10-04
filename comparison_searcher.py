#!/usr/bin/env python

#This python script reads in a list of desireable pairwise comparisons as a
#filtering step. The program then opens and searches a cuffdiff output file to
#identify transcripts/genes/cds files that have matching comparisons that were
#previously inputted in the first step. Usefully for shrinking a dataset to a
#size that is easier to manage in identifying any significant differences in
#expression.

#Files for reading including gene_exp.diff, isoform_exp.diff or cds_exp.diff or
#outouts from gene_expression_searcher.py

#Note all search comparison queries terms must be in same format as
#displayed in the cuffdiff output files.

#Code by Chris Gottschalk - 8/24/17 - email questions: gottsc33@msu.edu

import sys

#setting additional filtering for select pairwise comparisons
group = []
condition_details = raw_input("\n Please type in the $PATH to your list file of comparisons; enter each pairwise comparison must be in the same format as the -label option in cuffdiff seperated by a tab (ex. TRT1    CNT1) Press 'Enter' when done:")
in_compare = open(condition_details, "r")
with in_compare as comparison_list:
    for comparison in comparison_list:
        group.append(comparison.strip('\n'))
in_compare.close()
print "The following is the list you entered as pairwise comparisons you'd like to retrieve"
print group

gene_expression_searcher = raw_input("\n Please type in the $PATH to gene_expression_searcher or cuffdiff output Press 'Enter' when done:")
output1 = open(gene_expression_searcher, "r")
matches2 = []
for line in output1:
    matches_line = line.split()
    for comparison in group:
        comparison.strip('\t')
        if comparison in matches_line:
            matches2.append(line)
            break
output1.close()

output = open('comparison_searcher_filter.out', 'w+')
for line in matches2:
    output.write(line)
output.close()

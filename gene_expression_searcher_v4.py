#!/usr/bin/env python

#Version 4.0 - additional filtering for comparisons

#This python script reads in a list of transcript sequence identifiers that are
#of interest to you (a 1 column text file with a unique ID number one each line
#example AT1G05489 or MD17G099443). Usefully for shrinking a dataset to a size
#that is easier to manage in identifying any significant differences in
#expression.

#Version 2.1 upgraded a single line of code that could potential create slowed processing speeds due to a large input file being left open
#Version 3.0 integrated a iter count to keep track of the filtering and processing of the cuffdiff file
#Version 3.1 updated file upload instructions/code and prompts (retrieved from version 2.2 - unrealeased)
#Version 4.0 included secondary filtering for significant hits generating a secondary outfile

#Files for reading including gene_exp.diff, isoform_exp.diff, or cds_exp.diff.

#Note all search queries terms must be in same format as
#displayed in the cuffdiff output files.

#Code by Chris Gottschalk - 8/24/17 - email questions: gottsc33@msu.edu

import sys
import os
import time
import tqdm

#loading in target gene file

Welcome = """ 
&@                                                                @.                               
.@                                                                @                                
 @*   .......................................................    (@                                
 %@                                                              @#                                
  @@                                                            @@                                 
   @@      .............................................       @@                                  
    @@                                                        @@                                   
     .@(                                                    #@                                     
       (@&       .................................        @@*                                      
         .@@                                            @@                                         
            &@@                                      @@%                                           
               .@@@        .............         @@@                                               
                   /@@@                     .@@@*                                                  
                        %@@@            @@@%                                                       
                            ,@@@.  ,@@@                                                            
                               .@@@@.                                                              
                           .@@@      @@@                                                           
                        &@@.            %@@@@@@@@@@@@@@@                                           
                     @@@           *@@@*,,,,,,,,,,,,,,,,,&@@@                                      
                  @@(           @@@,,,,,,,*&@@@@@@@@@(,,,,,,,/@@                                   
               @@%            @@,,,,,*@@@(.,,,,,,,,,,,,@@@&,,,,,&@&                                
            .@@             @@,,,,(@@*,,,,,,,,,,,,,,,,,,,,,@@@,,,,/@#                              
          #@(     ........@@,,,,@@,,,,,,,,,,,,,,,,,,,,,,,,,,,,@@*,,,&@                             
        %@(              @&,,,&@*,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@@,,,,@&                           
      (@/              .@*,,,@@,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@(,,,@@                          
     @@    ............@*,,,@/,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@%,,,@%                         
    @#                @&,,,@(,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@(,,,@                         
  ,@                 *@,,,@@,,,,,,,,,,,,,,,,,,,,,,,,,@@,,,,,,,,,,,,,*@,,,@@                        
  @.  ,,,,,,,,,,,,,,,@@,,,@,,,,,,,,,,,,,,,,,,,,,,,@,&@,,,,,,,,,,,,,,,@@,,,@                        
 @&                  @/,,#@,,,,,,,,,,,,,,,,,,,,,,@,/@,,,,,,,,,,,,,,,,%@,,,@#                       
.@                   @*,,@@,,,,,,,,,,,,,,,,,,,,,@%,@,,,,,,,,,,,,,,,,,*@,,,@&                       
&@ ..................@*,,&@,,,,,,,,,,,,,,,,,,,,@@,@/,,,,,,,,,,,,,,,,,(@,,,@%                       
&@                   @&,,*@.,,,,,,,,,,,,,,,,,,@@,@%,,,,,,,,,,,,,,,,,,@@,,,@.                       
.@                   %@,,,@@,,,,,,,,,,,,,,,,,@@,@@,,,,,,,,,,,,,,,,,,.@*,,#@                        
 @/   ................@*,,*@,,,,,,,,,,,,,,,,*@,,%,,,,,,,,,,,,,,,,,,,@@,,,@(                        
 %@                   ,@,,,(@,,,,,,,,,,,,,,,@,,,,,,,,,,,,,,,,,,,,,,@@,,,&@                         
  @@                   &@,,,(@,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@@,,,%@                          
   @@      .............&@,,,,@@,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,@&,,,@@                           
    &@                    @&,,,*@@,,,,,,,,,,,,,,,,,,,,,,,,,,,,*@@,,,,@@                            
      @@                   @@,,,,*@@,,,,,,,,,,,,,,,,,,,,,,,,&@&,,,,@@                              
       ,@@       ,,,,,,,,,,,,@@*,,,,#@@%,,,,,,,,,,,,,,,,/@@@,,,,,@@.                               
          @@/                  (@@,,,,,,*@@@@@@@@@@@@@@(,,,,,,%@@                                  
            .@@#                  &@@#,,,,,,,,,,,,,,,,,,,,*@@@ @/                                  
                @@@        ...........#@@@@&/,,,,,,*#@@@@@      @@                                 
                    @@@%                    &@@&                *@@*,,*@@*                         
                        *@@@,          /@@@*                   @@,,,,,,,,@@                        
                             @@@%  &@@@                        @*,,,,,*@@(&@                       
                              .@@@@@@                          @&,,(@@*,,,,,@%                     
                           @@@.      .@@@                       @@@,,,,,,,,,,@@                    
                       @@@.              .@@@                    &@,,,,,,,@@@,@@                   
                    @@@                      @@@                  .@,,,@@@,,,,,*@*                 
                 @@(       .............        #@@                 @@@,,,,,,,,,,@@                
              @@(                                  #@@               @@,,,,,,,,,,,@@               
           ,@@                                        @@.             /@*,,,,,,,,,,/@              
         %@%     .................................      &@#             @(,,,,,,,,,,,@%            
       @@.                                                ,@&            @@,,,,,,,,,,,@@           
      @(                                                    #@            &@,,,,,,,,,,,#@          
    @@     .............................................      @@           .@(,,,,,,,,,,,@*        
   @@                                                          @@            @@,,,,,,,,,,,@@       
  @&                                                            @@            @@,,,,,,,,,,,@@      
 @@   ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,    @&            /@*,,,,,,,,,,/@,    
 @                                                               .@              @@,,,,,,,,,,,@&   
%@                                                                @               @@,,,,,,,,,,,@@  
 #                                                                %                ,@*,,,,,,,,,,/@.
                                                                                     @%,,,,,,,,,,,@
                                                                                      @@,,,,,,,,,,,
                                                                                       (@,,,,,,,,,,
                                                                                         @(,,,,,,,@
                                                                                          @@@&%@@@
                                                                                          
                                                                                           
Welcome to the interactive gene expression searcher - version 4.1

This software can be used to quickly parse cufflinks outputs to generate files that contain transcripts of interest.

NOTE: All paths for input files must be manually entered, unix autofill is not compatible.
"""
print(Welcome)
in_name = input("\n Please type in the $PATH to your target transcript sequence ID file and press 'Enter': ")
assert os.path.exists(in_name), "I did not find the file at, "+str(in_name)
in_gene = open(in_name, "r")
print("\n Your input file has been found and is loading...")
queries = []
with in_gene as query_list:
    for query in query_list:
        queries.append(query.strip('\n'))
in_gene.close()
#print queries #optional print list of query transcript sequence ID numbers

print("\n Input sequence IDs indexed...")

#reading in cuffdiff output

matches = []
cuffdiff_in = input("\n Please type in the $PATH to your cuffdiff output file and press 'Enter':")
assert os.path.exists(cuffdiff_in), "I did not find the file at, "+str(cuffdiff_in)
exp_data = open(cuffdiff_in, "r")
print("\n Your cuffdiff file has been found and is loaded, search starting...")

from tqdm import tqdm
for line in tqdm(exp_data): #first progress bar integration
    current_line = line.split()
    for query in queries: 
        if query in current_line:
            matches.append(line)
            #print matches
            break

#writing an output file that only contains information as origionally written by
#cuffdiff it will only contain lines with the unique transcript sequence ID
output = open('gene_expression_searcher.temp', 'w+')
print("\n Writing temp file of matches...")
for line in matches:
    output.write(line)
output.close()
exp_data.close()

#Setting additional filtering for select pairwise comparisons
comparison_matches = []
#User must define the comparison list below are the two need for the Gala & HC Transcriptomes
comparison_list = input("\n Please input all comparisons your would like to retrieve as written in cuffdiff ouput. example input format (include '' and regex tabs) - 'Control\tTreatment1, 'Control\tTreatment2' - press 'Enter when complete': ")
print("\n Thank you for entering comparisons, reading in temp file now...")
temp_in = open('gene_expression_searcher.temp', 'r')

with temp_in as non_compare:
    hits = non_compare.readlines()

print("\n Identifying target transcripts/genes and coorisponding desired comparisons...")

from tqdm import tqdm
for line in tqdm(hits): #integration of progress bar using tqdm
    for comparison in comparison_list:
        if comparison in line:
            comparison_matches.append(line)
            break
output2 = open('gene_expression_filter.txt', 'w+')
for line in comparison_matches:
    output2.write(line)  
output2.close()             
temp_in.close()
   
#Retrieve the significant DEgenes from filtered list to generated into a separate output
print("\n Identifying significantly DEgenes from filtered list...")  
sig_matches = []
significant = ['yes',]
from tqdm import tqdm
for line in tqdm(comparison_matches): #integration of progress bar using tqdm
    for phrase in significant:
        if phrase in line:
            sig_matches.append(line)
            break

print("\n Writing output file and removing temp files...")
 
	
sig_output = open('significant_filter.txt', 'w+')
for line in sig_matches:
    sig_output.write(line) 	
sig_output.close()	


#remove the temp file
os.remove('gene_expression_searcher.temp')
print("\n Thank you for choosing Gene Expression Searcher for all your cufflinks ouput needs! Happy transcriptome analyzing!")


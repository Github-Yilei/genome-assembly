#/usr/bin/python
#_*_coding: utf-8 _*_

"""
@ Description: the gaps of a genomes was filled with "N". 
@ Usage: python genome_gap_calculator.py input.file >output.file
@ Author: YiLei
@ Date: 2020-10-14
"""

import re
import sys
import argparse

def genome_gap_calculator(args):
    summary = 0
    i = 0
		results = "index\tchr_id\tstart\tend gap_length\tsummary\n"
    with open (args.input_path, "r") as genome_fa:
        for line in genome_fa:
						whole_line = ''
            if line.startswith(">"):
                chr_id = line.split()[0].replace(">", "")
            else:
                whole_line = line.replace("\n", "")
  
            for match in re.finditer("N+", whole_line):
                i += 1
                gap_length = match.end() - match.start()
                summary += gap_length
                results += str(i) + "\t" + chr_id + "\t" + str(match.start()) + "\t" + str(match.end()) + "\t" + str(gap_length) + "\t" + str(summary) + "\n"
            
    with open(args.output, 'w') as out_file:
        out_file.write(results )
        
  
if __name__=='__main__':
    parser = argparse.ArgumentParser(description =  "counting the length of gap and building a ploting-ready data frame for R")
    parser.add_argument('--input_path',
	                    dest = 'input_path',
						help = 'The absolute path of the genome fasta')
	
		parser.add_argument('--output', 
                        dest = "output", 
                        help = 'The result file in tsv format')
    
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    
    genome_gap_calculator(args)

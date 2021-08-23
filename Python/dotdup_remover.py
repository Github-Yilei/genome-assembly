#/usr/bin/python
#-*-coding: utf-8 -*-

"""
@ Description: The out put of gffread usually have duplicated protein sequence of one gene and  dot at the end of sequence those can not be recognisd by MCXscan or interproscan.
@ Usage: python dotdup_remover.py --input protein.fa --output_file' cleaned_protein.fa
@ Author: YiLei 2021.07.08
"""

import sys
import argparse

def dot_remover(args):
	gene_dict = {}
	with open(args.input, 'r') as lines:
		for line in lines:
			if line.startswith('>'):
				temp_line = ""
				gene_id = line.strip('>\n').split('.')[0]
				flag = gene_id
				try:
					gene_dict[gene_id]
				except KeyError:
					gene_dict[gene_id] = [line]
				else:
					gene_dict[gene_id].append(line)
			else:
				line = line.strip("\n")
				gene_dict[flag][-1] += line
    # select the longest sequence
	longest_seq = ""
	for k,v in gene_dict.items():
		if len(v) == 1:
			longest_seq += gene_dict[k][0].strip('.') + "\n"
		else:
			trans_max = ''
			for trans in gene_dict[k]:
				a = len(list(trans))
				b = len(list(trans_max))
				if a > b:
					trans_max = trans
				longest_seq += trans_max

	with open(args.output, 'w') as longest:
		longest.write(longest_seq)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This a protein sequence cleaner.')
    parser.add_argument('--input', dest = 'input', help = 'the protein sequence')
    parser.add_argument('--output_file', dest = 'output', help = 'the filtered contigs file')
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    dot_remover(args)

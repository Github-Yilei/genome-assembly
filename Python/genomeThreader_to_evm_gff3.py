#!/usr/bin/python
#-*- coding: utf-8 _*_

"""
@ Description: Preparing the ID-fields of the exons from genomeThreader which was necessary for EVidenceModeler to link up the alignment chains.
@ Usage: python genomeThreader_to_evm_gff3.py --input_file --output
@ Author: Yilei
@ Date:2020-01-05
"""

import re
import sys
import argparse

def genomeThreader_to_evm_gff3(args):
	exon = 0
	with open(args.input_file, "r") as gff:
		for line in gff:
			if not line.startswith("\n") and not line.startswith("#"):
				records = line.split("\t")
				records[1] = "gth"
			if re.search(r"\tgene\t", line):
				exon = 0
				tmp = records[8].replace('ID=', '')
				gene_id = re.sub(r';Target=.*', '', tmp).strip()
				records[8] = "ID={};Name={}".format(gene_id, gene_id)
			elif re.search(r"\texon\t", line):
				exon += 1
				exon_id  = gene_id + "." + str(exon)
				records[8] = "ID={};Parent={};Name={}".format(exon_id, gene_id, exon_id)
			else:
				continue
			with open(args.output, "a") as new_gff:
				new_gff.write("\t".join(records) +'\n')

if __name__=='__main__':
	parser = argparse.ArgumentParser(description = 'Rename the IDs of the gff3 file from EVidenceModeler')
	parser.add_argument('--input_file',
						dest = 'input_file',
						help = 'gff3 file from genomeThreader')
	parser.add_argument('--output', 
						dest = "output", 
						help = 'The result file')
											 
	if len(sys.argv) <= 1:
		parser.print_help()
		sys.exit()
	args = parser.parse_args()
    
	genomeThreader_to_evm_gff3(args)

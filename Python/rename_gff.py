#!/usr/bin/python
#-*- coding: utf-8 _*_

"""
@ Description: rename the IDs of the gff3 file from EVidenceModeler.
@ Usage: python rename_gff.py --input_file --prefix --output
@ Author: Yilei
@ Date:2020-01-13
"""

import re
import sys
import argparse

def rename_gff(args):
	count = 0
	mRNA  = 0
	cds   = 0
	exon  = 0
	prefix = args.prefix
	with open(args.input_file, "r") as gff:
		for line in gff:
			if not line.startswith("\n"):
				records = line.split("\t")
				records[1] = "."
			if re.search(r"\tgene\t", line):
				count = count + 10
				mRNA  = 0
				gene_id = prefix + str(count).zfill(6)
				records[8] = "ID={}".format(gene_id)
			elif re.search(r"\tmRNA\t", line):
				cds   = 0
				exon  = 0
				mRNA  = mRNA + 1
				mRNA_id  = gene_id + "." + str(mRNA)
				records[8] = "ID={};Parent={}".format(mRNA_id, gene_id)
			elif re.search(r"\texon\t", line):
				exon = exon + 1
				exon_id  = mRNA_id + "_exon_" + str(exon)
				records[8] = "ID={};Parent={}".format(exon_id, mRNA_id)
			elif re.search(r"\tCDS\t", line):
				cds = cds + 1
				cds_id  = mRNA_id + "_cds_" + str(cds)
				records[8] = "ID={};Parent={}".format(cds_id, mRNA_id)
			else:
				continue
			with open(args.output, "a") as new_gff:
				new_gff.write("\t".join(records) +'\n')

if __name__=='__main__':
	parser = argparse.ArgumentParser(description = 'Rename the IDs of the gff3 file from EVidenceModeler')
	parser.add_argument('--input_file',
						dest = 'input_file',
						help = 'gff3 file')
	parser.add_argument('--prefix',
						dest = 'prefix',
						help = 'the prefix of gene')
	parser.add_argument('--output', 
						dest = "output", 
						help = 'The result file')
											 
	if len(sys.argv) <= 1:
		parser.print_help()
		sys.exit()
	args = parser.parse_args()
    
	rename_gff(args)

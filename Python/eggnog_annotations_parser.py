#!/usr/bin/python
#-*- coding: utf-8 _*_

"""
@ Description: parsing the annotations file from eggnog-mapper and buidding a file for enrich analysis.
@ Usage: python eggnog_annotations_parser.py --input_file --key --output
@ Author: Yilei
@ Date:2020-02-06
"""

import sys
import argparse

def eggnog_annotations_parser(args):
	reslut = list()
	if args.key == 'GO':
		index = int(10) - 1
	elif args.key == 'KEGG':
		index = int(12) - 1
	elif args.key == 'info':
		index = int(9) - 1
	with open(args.input_file, "r") as tsv:
		for line in tsv:
			if line.startswith("#"):
				continue
			else:
				records = line.split("\t")
				geneid = records[0]
				tmp = records[index]
				tmp_split = tmp.split(',')
				for i in tmp_split:
					if i == '-':
						continue
					else:
						combined = geneid + "\t" + i + "\n"
						reslut += combined
		dup_reslut = list(set(reslut))
		sep = '\n'
		with open(args.output, "w") as term_file:
			#term_file.write(reslut)
			term_file.write(sep.join(dup_reslut))
		

if __name__=='__main__':
	parser = argparse.ArgumentParser(description = 'parsing the annotations file from eggnog-mapper and buidding a file')
	parser.add_argument('--input_file',
						dest = 'input_file',
						help = 'annotations file from eggnog-mapper')
						
	parser.add_argument('--key',
						dest = 'key',
						help = 'key words, can be GO or KEGG')					
	parser.add_argument('--output', 
						dest = "output", 
						help = 'The result file, plase check it by using "grep -o 'GO' HWB.emapper.annotations | wc -l"')
											 
	if len(sys.argv) <= 1:
		parser.print_help()
		sys.exit()
	args = parser.parse_args()
    
	eggnog_annotations_parser(args)

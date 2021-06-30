#!/usr/bin/python
#-*- coding: utf-8 _*_

"""
@ Description: parseing kegg json file and building K_to_ko.txt and K_to_ko for org.db making.
@ Url: https://www.kegg.jp/kegg-bin/get_htext?ko00001
@ Usage: python kegg_json_parser.py --input_file --output_path
@ Author: Yilei
@ Date:2020-02-08
"""
import sys
import argparse
import json

def kegg_json_parser(args):
	with open(args.input_file, 'r') as js:
		#ko_list = []
		K_to_ko = list()
		ko_to_Description = list()
    
		# convert json to dict
		map_dict = json.load(js)
		maps = map_dict['children']
		for layer_1 in maps:
			# https://www.genome.jp/kegg/pathway.html
			Pathway_Maps = layer_1['name'][6:]
        
			for layer_2 in layer_1['children']:
				# Pathway Identifiers
				Pathway_Identifiers = layer_2['name'][6:]
            
				for layer_3 in layer_2['children']:
					# this layer contain PATH information, [PATH:ko00010] 
					tmp = layer_3['name'][6:].split('[')
					ko_Description = tmp[0]
					try:
						tmp[1]
					except IndexError:
						continue
					else:
						ko_num =  tmp[1].strip(']').replace('PATH:', '').replace('BR:', '')
						combine_tmp = ko_num + "\t" + ko_Description
						ko_to_Description.append(combine_tmp)
                    
                    
					try:
						for layer_4 in layer_3['children']:
							# this is the last layer of kegg json, contain paired of gene and KO
							tmp = layer_4['name'].split('  ')
							K_num = tmp[0]
							K_Description =  tmp[1]
							K_ko = K_num + "\t" + ko_num
							K_to_ko.append(K_ko)

					except KeyError:
						continue
	sep = '\n'
	K2ko = args.outfile + '/K_to_ko.txt'
	with open(K2ko, 'w+') as K2ko:
		#term_file.write(reslut)
		K2ko.write(sep.join(K_to_ko))

	ko2description = args.outfile + '/ko_to_Description.txt'
	with open(ko2description, 'w+') as ko2description:
		ko2description.write(sep.join(ko_to_Description))
      
if __name__=='__main__':
	parser = argparse.ArgumentParser(description = 'parsing the annotations file from eggnog-mapper and buidding a file')
	
	parser.add_argument('--input_file',
						dest = 'input_file',
						help = 'annotations file from eggnog-mapper')	
	
	parser.add_argument('--output_path', 
						dest = "outfile", 
						help = 'The absolute path of the relsult file')
											 
	if len(sys.argv) <= 1:
		parser.print_help()
		sys.exit()
	args = parser.parse_args()
    
	kegg_json_parser(args) 

#/usr/bin/python
#_*_coding: utf-8 _*_

"""
@ Description: Parsing the json files of fastp and building a ploting-ready data frame for R.
@ Usage: python parse_fastp_json.py --input_path path --output ./result.txt
@ Author: YiLei
@ Date: 2021/06/16
"""

import sys
import pandas as pd
import argparse
from pathlib import Path

result_dict = {}
merged_result_dict = {}

def parse_fastp_json(args):
	p = Path(args,input_path)
	for child in p.iterdir():
		if child.endswith('json'):
			with open(child, 'r') as f:
				result_dict[child] = json.load(f)
				
				key = child.name.replace(".json", "")
				key1 = key + '+before_filtering'
				key2 = key + '+after_filtering'
				
				merged_result_dict[key1] = {'total_reads' : result_dict[i]["summary"]['before_filtering']['total_reads'],
                                            'total_bases' : result_dict[i]["summary"]['before_filtering']['total_bases'],
                                            'q20_rate' : result_dict[i]["summary"]['before_filtering']['q20_rate'],
                                            'q30_rate' : result_dict[i]["summary"]['before_filtering']['q30_rate']
                                            }
        merged_result_dict[key2] = {'total_reads' : result_dict[i]["summary"]['after_filtering']['total_reads'],
                                            'total_bases' : result_dict[i]["summary"]['after_filtering']['total_bases'],
                                            'q20_rate' : result_dict[i]["summary"]['after_filtering']['q20_rate'],
                                            'q30_rate' : result_dict[i]["summary"]['after_filtering']['q30_rate']
                                            }
							 
	df = pd.DataFrame(merge_result_dict).T
	df.to_csv(args.output, index = True, header = True)
	
	
if __name__=='__main__':
	parser = argparse.ArgumentParser(description =  "Parsing the josn files of fastp and building a ploting-ready data frame for R")
	
	parser.add_argument('--input_path',
	                    dest = 'input_path',
						          help = 'The absolute path of the fastp json file')
	
	parser.add_argument('--output', 
                        dest = "output", 
                        help = 'The result file in csv format')
	if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    
    parse_fastp_json(args)	

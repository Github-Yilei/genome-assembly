#/usr/bin/python
#code:utf-8
"""
@ Description: The GFF (General Feature Format) format consists of one line per feature,
@ each containing 9 columns of data, plus optional track definition lines.
@ Here, we try to build a a bed file for MCScanX .
@ Usage: python gfftoMCScanx.py  --input input_file --species species --output output_file
@ Date: 2020-10-15
"""
import sys
import argparse

def gfftoMCScanx(args):
	row_list = []
	with open(args.input, 'r') as anotation_lines:
		for line in anotation_lines:
			if line.startswith('#'):
				continue
			else:
				line_split = line.split()
				lines = dict()
				if line_split[0] == 'chrUn':
					continue
				elif elif line_split[2] != 'mRNA':
					continue
				else:
					whole_id = line_split[8]
					id_spl = whole_id.split(';')
					lines['ID'] = id_spl[0].replace('ID=', '')
					lines['Chr'] = line_split[0].lstrip('').replace('chr', chr(args.species)) # strip leading whitespace
					lines['start'] = int(line_split[3])
					lines['end'] = int(line_split[4])
				row_list.append(lines)
		data = pd.DataFrame(row_list)
		index = ["Chr", "ID", "start", "end"]
		result = data[index]
		result.to_csv(args.output, sep = "\t", header = False, index = False )
		
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This a protein sequence cleaner.')
    parser.add_argument('--input', dest = 'input', help = 'the protein sequence')
	parser.add_argument('--species', dest = 'species', help = 'a two-letter short name is used as prefix for the species')
    parser.add_argument('--output', dest = 'output', help = 'the filtered contigs file')
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    gfftoMCScanx(args)

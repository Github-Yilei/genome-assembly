#/usr/bin/python
#-*-coding: utf-8 -*-

"""
@ Description: This blastn parser dedacited in return the the largest pident, nident, or bitscore for plot.
@ Application: blastn-outfmt '6 qseqid sseqid pident nident qlen slen evalue bitscore'
@ Usage: python blast2plot.py --blast blastn.txt  --out myout.txt
@ Author: YiLei 2021.07.08
"""

import sys
import argparse

def blast2plot(args):
    pool = dict()
    with open(file, 'r') as lines:
        for line in lines:
            line_spl = line.split()
            flag = line_spl[0]
            try:
                pool[flag]
            except KeyError:
                pool[flag] = [line_spl]
            else:
                pool[flag].append(line_spl)

    my_result = ""
    for k,v in pool.items():
        if len(v) == 1:
            my_result += "\t".join(pool[k][0])
        else:
            tmp_list = list()
            tmp_max = 0
            for tmp in pool[k]:
                a = int(tmp[3])
                b = tmp_max
                if a > b:
                    tmp_max = a
                    tmp_list = tmp    
            my_result += "\t".join(pool[k][0])

    with open(output_file, 'w') as longest:
        longest.write(my_result)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'This a blastn parser.')
    parser.add_argument('--blast', dest = 'blast', help = 'the result of blastn')
    parser.add_argument('--out', dest = 'output_file', help = 'the filtered contigs file')
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    blast2plot(args) 

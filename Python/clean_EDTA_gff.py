"""
@ Description: after EDTA pip, the repetitive elements is in /short_name_genome.mod.EDTA.TElib.fa
@ Description: this script is to transform the short_name_genome.mod.EDTA.TElib.fa for further stats
@ Usage: python clean_TElib_fa.py input_file output_file
@ file: clean_TElib_fa.py
@ Time: 2020-10-18
"""

import sys
import re
input_file = sys.argv[1]
output_file = sys.argv[2]


##
clean_read = ''

with open(input_file, 'r') as TElib_gff:
    for line in TElib_gff:
        if line.startswith('#'):
            continue
        else:
            line_spl = line.split()
            start = line_spl[3]
            end = line_spl[4]
            attributes = line_spl[8]
            length = int(end) - int(start)
#            attributes_spl = attributes.split(';')
            classification = ''.join(re.findall(r"Classification=(.+?);", attributes))
#            classification = attributes_spl[2].replace('Classification=', '')
            classification_spl = classification.split('/')

            try:
                 classification_spl[1]
            except(IndexError, ValueError, ArithmeticError):
                 classification_spl.append("Missing")
            finally:
                LTRs =  classification_spl[0].replace('\n', '')
                LTR_RTS =  classification_spl[1]
                LTR_RTS = LTR_RTS.replace('\n', '')
            line_clean = str(length) + '\t' + LTRs + '\t' + LTR_RTS + '\t' + attributes + '\t' + end + '\n'
        # combine reslut
        clean_read += line_clean
        # export reads
with open(output_file, 'w') as output:
    output.write(clean_read)

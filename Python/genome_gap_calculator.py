"""
@ Description: the gaps of a genomes was filled with "N". 
@ Usage: python genome_gap_calculator.py input.file >output.file
@ Author: YiLei
@ Date: 2020-10-14
"""


import re
import sys
import_file = sys.argv[1]

summary = 0
lines = ""
i = 0

with open (import_file, "r") as genome_fa:
    for line in genome_fa:
        if line.startswith(">"):
            line = line.split()[0]
        else:
            line = line.replace("\n", "")
        lines += line

    for match in re.finditer("N+", lines):
        i += 1
        gap_length = match.end() - match.start()
        summary = summary + gap_length
        
        print(i, match.start(), match.end(), gap_length, summary)

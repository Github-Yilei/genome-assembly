"""
@ Description: this secripts take MCScanX_collinearity and MCScanX_gff as input, and return the circos collinearity data frame.
@ Usage: python make_circos_links.py MCScanX_collinearity MCScanX_gff out_put
"""

import pandas as pd
import sys

collinearity = sys.argv[1]
gff = sys.argv[2]
output_file = sys.argv[3]

collinearity_data = list()
with open(collinearity, 'r') as coll:
    for line in coll:
        index = {}
        if line.startswith('#'):
            continue
        else:
            line_spl = line.split()
            index["ID1"] = line_spl[2]
            index["ID2"] = line_spl[3]
            collinearity_data.append(index)

collinearity_df = pd.DataFrame(collinearity_data)


# set gff
gff_list = list()

with open(gff, 'r') as gff:
    for lines in gff:
        tmp_gff = {}
        line_spl = lines.split()
        tmp_gff['chr'] = line_spl[0]
        tmp_gff['ID'] = line_spl[1]
        tmp_gff['start'] = line_spl[2]
        tmp_gff['end'] = line_spl[3]
        gff_list.append(tmp_gff)
gff_df = pd.DataFrame(gff_list)

# reorder for rename
index = ['chr', 'ID', 'start', 'end']
gff_df = gff_df[index]

# combine species1
col_species1 = ['chr1', 'ID1', 'start1', 'end1']
gff_df.columns = col_species1
tmp_merge = pd.merge(collinearity_df, gff_df, how='left', on= 'ID1')
#print(tmp_merge)

# combine species2
col_species2 = ['chr2', 'ID2', 'start2', 'end2']
gff_df.columns = col_species2
result_merge = pd.merge(tmp_merge, gff_df, how='left', on= 'ID2')

# set data for circos
index = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']
reslut = result_merge[index]

reslut.to_csv(output_file, sep = "\t", header = False, index = False )

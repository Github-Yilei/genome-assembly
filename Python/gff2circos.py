
"""
@ Description: The GFF (General Feature Format) format consists of one line per feature,
@ Description: each containing 9 columns of data, plus optional track definition lines.
@ Description: this sript will extract information from .gff3 with key words: CDS, intron, exon, gene, etc.
@ then, build a circos windows density data frame with chromosome lenght file
@ Usage: python circos_windows_anotation_gff.py gff_file chr_length keywords windows_step circos_outfile
@ file: circos_windows_anotation_gff.py
@ Date: 2020-10-15
@ Check the output by using awk: awk '{a+=$4}END{print a}' circos_outfile the output should be the same as data row number.
"""

import pandas as pd
import sys



input_file = sys.argv[1]
input_chr_length = sys.argv[2]
keywords = sys.argv[3]
windows_step = sys.argv[4]
output_file = sys.argv[5]


row_list = []
windows = int(windows_step)

with open(input_file, 'r') as anotation_lines:
    for line in anotation_lines:
        if line.startswith('#'):
            continue
        else:
            line_split = line.split('\t')
            lines = dict()
            if line_split[2] != str(keywords):
                continue
            else:
                lines['Chr'] = line_split[0].lstrip('') # strip leading whitespace
            #lines['source'] = line_split[1]
            #lines['feature'] = line_split[2]
                lines['start'] = int(line_split[3])
                lines['end'] = int(line_split[4])
            #line['score'] = line_split[5]
            #lines['length'] = int(line_split[4])-int(line_split[3])
            #lines['strand'] = line_split[6]
            #lines['frame'] = line_split[7]
            row_list.append(lines)

##############################
# extract start end
start = []
end = []

for line in row_list:
    start.append(line['start'])
    end.append(line['end'])

# set index
max_range = max(end)/windows
if isinstance(max_range, float) == True:
    max_range = int(max_range) + 1

index = range(0, max_range + 1)

# set window_start window_end
window_start = []
window_end = []

for i in end:
    for j in index:
        if i >= (j)*windows and i < (j+1)*windows:
            window_start.append((j)*windows)
            window_end.append((j+1)*windows)

# combine result
data = pd.DataFrame(row_list)
data["window_start"] = window_start
data["window_end"] = window_end

# print(data)

# merge chr length to data
chr_list = []
with open (input_chr_length) as chr_data:
    for lines in chr_data:
        line = lines.replace('\n', '')
        temp = line.split('\t')
        chr_list.append(temp )

chr_df = pd.DataFrame(chr_list, columns=['Chr','chr_length'])

data_merged = pd.merge(data, chr_df, how='left', on= 'Chr')

print("notice：\nthe number of rows are:", len(data_merged), "\n")
print("please check the output with: awk '{a+=$4}END{print a}' out_file\n")
print("the output should be the same as the data rows\n\n")
# build index
# index = ['Chr', 'feature', 'start', 'end', 'score', 'frame', 'strand']
# data = data[index]
# print(data)
def compare_return(window_start, window_end, chr_length):
    chr_length = int(chr_length)
    if window_start < chr_length and chr_length <= window_end:
        new_window_end = chr_length
    else:
        new_window_end = window_end
    return(new_window_end)


data_merged['new_end'] = data_merged.apply(lambda x: compare_return(x["window_start"], x["window_end"], x["chr_length"]), axis = 1)

# print(data_merged)
result = data_merged.loc[ :,["Chr", "window_start", "new_end"]]

# density statistics
# colnames of the results is the same as the counted colums
# if useing size() the colname is 0
df = result.groupby(['Chr', 'window_start', "new_end"]).size() #.count()
df = df.reset_index()

#df["end"] = result["new_window_end"] for count()

index = ['Chr', 'window_start', 'new_end', 0]
df = df[index]

df.to_csv(output_file, sep = "\t", header = False, index = False )

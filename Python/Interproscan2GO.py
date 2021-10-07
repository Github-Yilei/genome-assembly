import sys
import re

input_file = sys.argv[1]
out = sys.argv[2]

GO_list = list()
with open(input_file) as lines:
    data = dict()
    for line in lines:
        line_re = re.search(r".+GO:.+", line)
        if line_re == None:
            continue
        else:
            line_match = line_re.group()

            line_match_spl = line_match.split("\t")
            geneid = line_match_spl[0]
            GOs = line_match_spl[13]
            GOs_spl = GOs.split("|")
            for i in GOs_spl:
                 combined = geneid + "\t" + i + "\t" + "IEA"
                 GO_list.append(combined)
#            data.setdefault(index, 'GO')
dup_GO_list = list(set(GO_list))
dup_GO_list.insert(0, "GID\tGO\tEVIDENCE")
sep = "\n"

with open(out, 'w+') as GID2GO:
    GID2GO.write(sep.join(dup_GO_list))

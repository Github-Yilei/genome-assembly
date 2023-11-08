import sys
import pathlib
# file_path：the ancestor path of EDTA result
file_path = sys.argv[1]

sample_idx = {}
p = pathlib.Path(file_path)
for child in p.iterdir():
    if child.is_dir():
        for sub_child in child.iterdir():
            if 'mod.EDTA.TEanno.sum' in sub_child.name:
                sample_idx[sub_child.parent.name] = sub_child


for sample in sample_idx.keys():
    with sample_idx[sample].open('r') as lines:
        i = 0
        for line in lines:
            line_spl = line.strip().split()
            i += 1
            # LTR
            if i == 8:
                print(sample + "\t" + "LTR\tCopia\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 9:
                print(sample + "\t" + "LTR\tGypsy\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 10:
                print(sample + "\t" + "LTR\tunknown\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # TIR
            if i == 12:
                print(sample + "\t" + "TIR\tCACTA\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 13:
                print(sample + "\t" + "TIR\tMutator\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 14:
                print(sample + "\t" + "TIR\tPIF_Harbinger\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 15:
                print(sample + "\t" + "TIR\tTc1_Mariner\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 16:
                print(sample + "\t" + "TIR\thAT\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # nonLTR
            if i == 18:
                print(sample + "\t" + "nonLTR\tLINE_element\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if i == 19:
                print(sample + "\t" + "nonLTR\tunknown\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # nonTIR
            if i == 21:
                print(sample + "\t" + "nonTIR\thelitron\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # repeat_region
            if i == 22:
                print(sample + "\t" + "repeat_region\trepeat_region\t" + "\t".join(str(e) for e in  line_spl[-3:]))

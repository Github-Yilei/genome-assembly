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
            if "Copia" in line_spl:
                print(sample + "\t" + "ClassI\tLTR\tCopia\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "Gypsy" in line_spl:
                print(sample + "\t" + "ClassI\tLTR\tGypsy\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "unknown" in line_spl and i < 15:
                print(sample + "\t" + "ClassI\tLTR\tunknown\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # TIR
            if "CACTA" in line_spl:
                print(sample + "\t" + "ClassII\tTIR\tCACTA\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "Mutator" in line_spl:
                print(sample + "\t" + "ClassII\tTIR\tMutator\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "PIF_Harbinger" in line_spl:
                print(sample + "\t" + "ClassII\tTIR\tPIF_Harbinger\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "Tc1_Mariner" in line_spl:
                print(sample + "\t" + "ClassII\tTIR\tTc1_Mariner\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "hAT" in line_spl:
                print(sample + "\t" + "ClassII\tTIR\thAT\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            #if "polinton" in line_spl:    
            #    print(sample + "\t" + "ClassII\tTIR\tpolinton\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # nonLTR
            if "LINE_element" in line_spl:
                print(sample + "\t" + "ClassI\tnonLTR\tLINE_element\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            if "unknown" in line_spl and i > 15:
                print(sample + "\t" + "ClassI\tnonLTR\tunknown\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # nonTIR
            if "thelitron" in line_spl:
                print(sample + "\t" + "ClassII\tnonTIR\thelitron\t" + "\t".join(str(e) for e in  line_spl[-3:]))
            # repeat_region
            if "repeat_region" in line_spl:
                print(sample + "\t" + "Repeat\trepeat_region\trepeat_region\t" + "\t".join(str(e) for e in  line_spl[-3:]))



import sys
import argparse

def ChromWindow(fai, windows):
   chrom_windows = list()
   with open(fai, 'r') as lines:
        for line in lines:
            line = line.strip()
            line_spl = line.split()
            chrom = line_spl[0]
            end = int(line_spl[1])
            windows = int(windows)
            max_index = int(end/windows)

            # buliding windows
            steps = []
            for i in range(0, end, windows):
                steps.append(i)

            if max_index * windows == end:
                start_temp = steps

                end_temp = steps[1:]
                end_temp.append(end)

            elif max_index * windows < end:
                start_temp = steps
                start_temp.append((max_index + 1) * windows)

                end_temp = start_temp[1: ]
                end_temp.append(end)

            for i in range(0, len(end_temp)):
                chrom_windows.append([chrom, start_temp[i], end_temp[i]])

        return chrom_windows



def HapvsHet(args):
    windows_position = ChromWindow(args.fai, args.windows)
    sample_idx = []

    hapvshet_mat = dict()
    with open(args.input_tsv, 'r') as mat:
        for line in mat:
            line = line.strip()
            if line.startswith('#CHROM'):
                sample_idx = line.split("\t")

            else:
                line_spl = line.split("\t")
                chrom = line_spl[0]
                pos = line_spl[1]

                for i in range(2, len(line_spl)):
                    type_of_snp = line_spl[i]
                    sample_name =  sample_idx[i]

                    for pos_list in windows_position:
                        if chrom == pos_list[0] and int(pos) >= int(pos_list[1]) and int(pos) <= int(pos_list[2]):
                            key = pos_list[0] + '_' + str(pos_list[1]) + '_' + str(pos_list[2]) + '_' + sample_name
                            try:
                                hapvshet_mat[key]
                            except KeyError:
                            #   ref, het, qry, mis, total
                                hapvshet_mat[key] = [0, 0, 0, 0, 0, key]
                            else:
                                hapvshet_mat[key][4] = hapvshet_mat[key][4] + 1
                                if type_of_snp == 'ref':
                                    hapvshet_mat[key][0] += 1
                                elif type_of_snp == 'het':
                                    hapvshet_mat[key][1] +=  1
                                elif type_of_snp == 'qry':
                                    hapvshet_mat[key][2] += 1
                                elif type_of_snp == 'mis':
                                    hapvshet_mat[key][3] +=  1

                   # print(hapvshet_mat)
    with open(args.output, 'w') as out_put:
        out_put.write("hap" + "\t" + "het" + "\t" + "qry" + "\t" + "mis" +  "\t" + "total" + "\t" + "sample" + "\n")
        for value in hapvshet_mat.values():
            out_put.write(str(value[0]) + "\t" + str(value[1]) + "\t" + str(value[2]) + "\t" +
                          str(value[3]) + "\t" + str(value[4]) + "\t" + value[5] + "\n")

if __name__=='__main__':
        parser = argparse.ArgumentParser(description = 'parsing the tsv from vcf2matrix.py')
        parser.add_argument('--fai',
                            dest = 'fai',
                            help = 'fai file of samtools index')
        parser.add_argument('--windows',
                            dest = 'windows',
                            help = 'windows step of chromsome')

        parser.add_argument('--input',
                            dest = 'input_tsv',
                            help = 'tsv file')

        parser.add_argument('--output',
                            dest = "output",
                            help = 'The result file')

        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        args = parser.parse_args()

        HapvsHet(args)

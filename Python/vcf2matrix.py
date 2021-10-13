import sys
import argparse

def SampleList(sample_file):
    sample_list = list()
    with open(sample_file, 'r') as samples:
        for line in samples:
            id = line.strip()
            sample_list.append(id)

    return sample_list
    
    
def vcf2matrix(args):
    sample_list = SampleList(args.sample_file)
    snp_list = list()
    sample_idx = list()
    with open(args.input_vcf, 'r') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue
            elif line.startswith('#C'):
                line = line.strip()
                line_spl = line.split("\t")
                sample_id = line_spl[0:2]
                for i in sample_list:
                    sample_idx.append(line_spl.index(i))
                    sample_id.append(line_spl[line_spl.index(i)])
                sample_id.append('ref_count')
                snp_list.append(sample_id)
            else:
                temp = list()
                ref_count = 0
                line_spl = line.split("\t")
                temp += line_spl[0:2]
                for i in sample_idx:
                    GT = line_spl[i]
                    if '0/0' in GT:
                        temp.append('ref')
                        ref_count += 1
                    elif '0/1' in GT:
                        temp.append('het')
                    elif '1/1' in GT:
                        temp.append('qry')
                    elif './.' in GT:
                        temp.append('mis')
                    else:
                        temp.append('ukn')
                temp.append(str(ref_count))
                snp_list.append(temp)

    with open(args.output, 'a') as f:
        for lst in snp_list:
            f.write("\t".join(lst) + "\n")

if __name__=='__main__':
        parser = argparse.ArgumentParser(description = 'parsing the vcf and building a matrix for identity-by-descent analysis')
        parser.add_argument('--samples', 
                            dest = 'sample_file', 
                            help = 'sample names file, one name per line')
                            
        parser.add_argument('--input', 
                            dest = 'input_vcf', 
                            help = 'vcf file')
                            
        parser.add_argument('--output', 
                            dest = "output", 
                            help = 'The result file')

        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        args = parser.parse_args()

        vcf2matrix(args)

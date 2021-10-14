import sys
import argparse
import re

def GenePosition(gff3):
    gene_position = list()
    with open(gff3, 'r') as gff:
        for line in gff:
            line = line.strip()
            line_spl = line.split()
            if line_spl[2] == 'gene':
                chrom = line_spl[0:1]
                start_end = line_spl[3:5]
                gene_len = int(line_spl[4]) - int(line_spl[3])
                name = re.sub(".*Name=", "", line_spl[8])
                gene_position.append(chrom + start_end + [name] + [str(gene_len)])

    return gene_position


def SnpPerGene(args):
    gene_position = GenePosition(args.gff3)
    conuter_dict = dict()
    with open(args.input_vcf, 'r') as vcf:
        for line in vcf:
        #    line = line.strip()
            if line.startswith('#'):
                continue
            else:
                line_spl = line.split("\t")
                chrom = line_spl[0]
                pos = line_spl[1]
                for pos_lis in gene_position:
                    if chrom == pos_lis[0] and int(pos) >= int(pos_lis[1]) and int(pos) <= int(pos_lis[2]):
                        try:
                            conuter_dict[pos_lis[3]]
                        except KeyError:
                            conuter_dict[pos_lis[3]] = 0
                        else:
                            conuter_dict[pos_lis[3]] += 1
                        break

    for pos_lis in gene_position:
        gene_id = pos_lis[3]
        gene_len = int(pos_lis[4])
        try:
            conuter_dict[gene_id]
        except KeyError:
            num_of_snp = 0
        else:
            num_of_snp = int(conuter_dict[gene_id])
        num_per_bp = "{}".format(num_of_snp/gene_len)
        with open(args.output, 'a') as f:
            f.write(gene_id + "\t" + num_per_bp + "\t" + str(num_of_snp) + "\t" + str(gene_len) + "\n")

if __name__=='__main__':
        parser = argparse.ArgumentParser(description = 'parsing the vcf and building a matrix with the number of SNPs per gene')
        parser.add_argument('--gff',
                            dest = 'gff3',
                            help = 'gff3 file')

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

        SnpPerGene(args)

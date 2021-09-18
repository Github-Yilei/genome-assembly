import sys
import argparse
from pathlib import Path

def multfa2OG(args):
        path = Path(args.path)
        gene_dict = {}
        flag = ''
        with open(args.multfa, 'r') as multfa:
                for line in multfa:
                        if line.startswith('>'):
                                gene_id = line.strip('>\n').split()[0]
                                flag = gene_id
                                try:
                                        gene_dict[gene_id]
                                except KeyError:
                                        gene_dict[gene_id] = ['>' + gene_id + '\n']
                                else:
                                        print("ERROR: Gene id is Duplicated, please check your mult-fasta file.")
                                        sys.exit()
                        else:
                                # adding sequence at the end of gene_id
                                gene_dict[flag][-1] += line

        with open(args.OGs, 'r') as lines:
                for line in lines:
                        temp_fasta = str()
                        line_spl = line.strip().split(':')
                        OG_id = line_spl[0] + '.fa'
                        gene_id = line_spl[1].split()
                        for key in gene_id:
                                try:
                                        temp_fasta += gene_dict[key][0]
                                except KeyError:
                                        continue
                        temp = path.joinpath(OG_id)
                        with open(temp, 'w+') as OG_fasta:
                                OG_fasta.write(temp_fasta)

if __name__=='__main__':
        parser = argparse.ArgumentParser(description = 'extracting fasta sequence from a mult-fasta file.')

        parser.add_argument('--OGs',
                                                dest = 'OGs',
                                                help = 'a file likes Orthogroups.txt with OGs and gene id.')

        parser.add_argument('--multfa',
                                                dest = 'multfa',
                                                help = 'a combined fasta file.')

        parser.add_argument('--path',
                                                dest = 'path',
                                                help = 'The absolute path of the results.')
        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        args = parser.parse_args()

        multfa2OG(args)

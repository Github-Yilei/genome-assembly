import subprocess
import sys

one2one = sys.argv[1]
bed = sys.argv[2]
singleton = sys.argv[3]

anchor_list = []
with open(one2one, 'r') as anchors:
        for line in anchors:
                if line.startswith('#'):
                        continue
                else:
                        line_spl = line.split()
                        anchor_list.append(line_spl[0])


gene_order_list = []
hap_bed = dict()
with open(bed, 'r') as bed:
        for line in bed:
                line_spl = line.split()
                gene_order_list.append(line_spl[3])
                key = line_spl[3]
                try:
                        hap_bed[key]
                except KeyError:
                        hap_bed[key] = line_spl

with open(singleton, 'r') as singletons:
        for gene in singletons:
                gene = gene.strip()
                gene_idx = gene_order_list.index(gene)

                # absense gene
                absense_gene = gene_order_list[gene_idx]
                absense_gene_chr =  hap_bed[absense_gene][0]

        # searching for neighbor anchor
                for i in range(1,51):
                        try:
                                temp = gene_order_list[gene_idx - i]
                        except IndexError:
                                front_gene = gene_order_list[gene_idx - i + 1]
                        else:
                                if hap_bed[temp][0] != absense_gene_chr:
                                        front_gene = gene_order_list[gene_idx - i + 1]
                                        break
                                elif temp in anchor_list:
                                        front_gene =  gene_order_list[gene_idx - i]
                                        break

                F_start = hap_bed[front_gene][0:2]
                F_gene = hap_bed[front_gene][3]

                for i in range(1,51):
                        try:
                                temp = gene_order_list[gene_idx + i]
                        except IndexError:
                                end_gene = gene_order_list[gene_idx + i -1]
                        else:
                                if hap_bed[temp][0] != absense_gene_chr:
                                        end_gene = gene_order_list[gene_idx + i -1]
                                        break
                                elif temp in anchor_list:
                                        end_gene = gene_order_list[gene_idx + i]
                                        break

                R_end = hap_bed[end_gene][2:3]
                R_gene = hap_bed[end_gene][3]

                # building new bed
                gene_pair = F_gene + ":" + R_gene
                bed_line = F_start + R_end
                bed_line.append(gene)
                bed_line.append(gene_pair)

                print("\t".join(bed_line))
                #write(hap1_hap2.singleton.bed)
# subprocess.run("~/miniconda3/bin/seqkit subseq --bed hap1_hap2.singleton.bed hap1_chrom.fa -o hap1_hap2.singleton.fa", shell=True)
# sed -i "s#hap.*\:\. ##" hap1_hap2.singleton.fa
# ~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/makeblastdb -in chrom.fa -parse_seqids -dbtype nucl -out chrom_db
# ~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/blastn -num_threads 10 -evalue 1e-10 -outfmt '6 qseqid sseqid pident nident qlen slen evalue bitscore'  -db hap2_db -query hap1_hap2.singleton.fa -out hap1_hap2.blast
#  awk -F"\t" '$3>=95  {print $1}' hap1_hap2.blast.score  | sort | uniq | wc -l
#

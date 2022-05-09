import sys
import re
input_file = sys.argv[1]
pep_out = sys.argv[2]
gff_out = sys.argv[3]

gff_list =  list()
uniq_id_list = list()
with open(input_file, 'r') as gff3:
    temp_pep = dict()
    for line in gff3:
        if line == '\n':
            continue
        elif line.startswith('# ORIGINAL') or line.startswith('# PASA_UPDATE'):
            temp_pep = dict()
        elif line.startswith('chr'):
            temp = dict()
            ID = ''
            line_spl = line.split("\t")
            line_spl[1] = 'PASA'
            chrom = line_spl[0]
            start = line_spl[3]
            end = line_spl[4]

            if line_spl[2] == 'gene':
                ID = 'a'
                temp_pep['chrom'] = chrom
                temp_pep['start'] = start
                temp_pep['end'] = end
                temp_pep['ID'] = 'g'
                temp_pep['GFF'] = ''
            elif line_spl[2] == 'mRNA':
                ID = 'b'
            elif line_spl[2] == 'five_prime_UTR':
                ID = 'c'
            elif line_spl[2] == 'exon':
                ID = 'd'
            elif line_spl[2] == 'CDS':
                ID = 'e'
            elif line_spl[2] == 'three_prime_UTR':
                ID = 'f'

            temp['chrom'] = chrom
            temp['start'] = start
            temp['end'] = end
            temp['ID'] = ID
            temp['GFF'] = '\t'.join(line_spl[:-1]) + "\tPlaceholders"
            # uniq lines
            if temp['GFF'] not in uniq_id_list:
                gff_list.append(temp)
                uniq_id_list.append(temp['GFF'])
        elif line.startswith('#PROT'):
            pep_spl = line.split("\t")
            pep_seq = pep_spl[1]
            # if this is the first prot
            if len(temp_pep['GFF']) == 0:
                temp_pep['GFF'] = 'represnt_pep' + "\t" + pep_seq
                gff_list.append(temp_pep)

            elif len(pep_seq) > len(temp_pep['GFF']):
                temp_pep['GFF'] = 'represnt_pep' + "\t" + pep_seq
                gff_list[-1] = temp_pep

sorted_list = sorted(gff_list, key=lambda k: (k['chrom'], int(k['start']), k['ID']))

count = 0
mRNA  = 0
cds   = 0
exon  = 0
UTR_5 = 0
UTR_3 = 0
prefix = "Cp"

for i in range(len(sorted_list)):
        temp = sorted_list[i]
        line = temp['GFF'].strip()
        # renamer
        if  line.startswith('chr'):
            records = line.split("\t")
            records[1] = "PASA"
        if re.search(r"\tgene\t", line):
            count = count + 10
            mRNA  = 0
            UTR_5 = 0
            UTR_3 = 0
            chr_num = records[0]
            #gene_id = chr_num + '_' + str(count).zfill(7)
            gene_id = prefix + '_' + chr_num + '_' + str(count).zfill(6)
            pep_id = ">" + gene_id + "\t" + "gene=" + gene_id + "\n"
            records[8] = "ID={};Name={}".format(gene_id, gene_id)

        elif line.startswith('represnt_pep'):
            pep_spl = line.split("\t")
            pep_seq = pep_spl[1].strip("*")
            pep_records = pep_id + pep_seq  + "\n"
            with open(pep_out, "a") as pep_file:
                pep_file.write(pep_records)
            pep_id = ''
            records = ''
        elif re.search(r"\tmRNA\t", line):
            cds   = 0
            exon  = 0
            mRNA  = mRNA + 1
            mRNA_id  = gene_id + "." + str(mRNA)
            records[8] = "ID={};Parent={};Name={}".format(mRNA_id, gene_id, mRNA_id)
        elif re.search(r"\texon\t", line):
            exon = exon + 1
            exon_id  = mRNA_id + "_exon_" + str(exon)
            records[8] = "ID={};Parent={};Name={}".format(exon_id, mRNA_id, exon_id)
        elif re.search(r"\tCDS\t", line):
            cds = cds + 1
            cds_id  = mRNA_id + "_cds_" + str(cds)
            records[8] = "ID={};Parent={};Name={}".format(cds_id, mRNA_id, cds_id)
        elif re.search(r"\tfive_prime_UTR\t", line):
            UTR_5 = UTR_5 + 1
            UTR_5_id = gene_id + ".UTR_5." + str(UTR_5)
            records[8] = "ID={};Parent={}".format(UTR_5_id, gene_id)
        elif re.search(r"\tthree_prime_UTR\t", line):
            UTR_3 = UTR_3 + 1
            UTR_3_id = gene_id + ".UTR_3." + str(UTR_3)
            records[8] = "ID={};Parent={}".format(UTR_3_id, gene_id)
        else:
            continue
        with open(gff_out, "a") as new_gff:
            new_gff.write("\t".join(records) +'\n')


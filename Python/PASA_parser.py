'''
PASA_parser was used to parse the gff3 file of Gene Structure Annotation and Analysis Using PASA and reslut into sorted gff3 and pep file.
'''
import sys
import re
input_file = sys.argv[1]
pep_out = sys.argv[2]
gff_out = sys.argv[3]

gff_dict =  dict()
uniq_id_list = list()

with open(input_file, 'r') as gff3:
    for line in gff3:
        if line == '\n':
            continue
        elif line.startswith('# ORIGINAL') or line.startswith('# PASA_UPDATE'):
            gff_key = ''
            temp_pep = ''
            temp = {}
        elif line.startswith('chr'):
            line_spl = line.split("\t")
            line_spl[1] = 'PASA'
            gff_value = '\t'.join(line_spl[:-1]) + "\tPlaceholders"

            if line_spl[2] == 'gene':
                chrom = line_spl[0]
                start = line_spl[3]
                end = line_spl[4]

                gff_key = chrom + '_' + start +'_' + end
                temp = {'chrom' : chrom, 'start' : start, 'end' : end}

                uniq_id_list.append(temp)
                gff_dict[gff_key] = [gff_value]
            else:
                gff_dict[gff_key].append(gff_value)

        elif line.startswith('#PROT'):
            pep_spl = line.split("\t")
            pep_seq = pep_spl[1]
            # if this is the first prot
            if len(temp_pep) == 0:
                temp_pep = 'represnt_pep' + "\t" + pep_seq
                gff_dict[gff_key].append(temp_pep)

            elif len(pep_seq) > len(temp_pep):
                temp_pep = 'represnt_pep' + "\t" + pep_seq
                gff_dict[gff_key][-1] = temp_pep
# sort keys
sorted_list = sorted(uniq_id_list, key=lambda k: (k['chrom'], int(k['start'])))

count = 0
mRNA  = 0
cds   = 0
exon  = 0
UTR_5 = 0
UTR_3 = 0
prefix = "Cp"

for i in range(len(sorted_list)):
    temp = sorted_list[i]
    gff_key = temp['chrom'] + '_' + temp['start'] + '_' + temp['end']
    temp_list = gff_dict[gff_key]
    for i in range(len(temp_list)):
        line = temp_list[i]
        records = line.split("\t")
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
            pep_records = pep_id + pep_seq
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
        # skip the position of pep seq
        if len(records ) == 9:
            with open(gff_out, "a") as new_gff:
                new_gff.write("\t".join(records) +'\n')

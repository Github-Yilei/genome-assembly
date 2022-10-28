import sys
import argparse
import subprocess

def runcmd(command):
    runs = subprocess.run(command, shell=True, encoding="utf-8")
    if runs.returncode == 0:
        print("success:",runs)
    else:
        print("error:",runs)

def perform_blast(args):
    seq_1 = args.i1
    seq_2 = args.i2
    seq_type = args.t
    if seq_type == "n":
        mkdb1 = "makeblastdb -dbtype nucl -input_type fasta -out seq_1.db -in "  + seq_1
        mkdb2 = "makeblastdb -dbtype nucl -input_type fasta -out seq_2.db -in "  + seq_2

        blast1 = "blastn -db seq_2.db -outfmt 6 -out seq_1.out -query "  + seq_1
        blast2 = "blastn -db seq_1.db -outfmt 6 -out seq_2.out -query "  + seq_2


    elif seq_type == "p":
        mkdb1 = "makeblastdb -dbtype prot -input_type fasta -out seq_1.db -in " + seq_1
        mkdb2 = "makeblastdb -dbtype prot -input_type fasta -out seq_2.db -in " + seq_2

        blast1 = "blastp -db seq_2.db -outfmt 6 -out seq_1.out -query "  + seq_1
        blast2 = "blastp -db seq_1.db -outfmt 6 -out seq_2.out -query "  + seq_2

    runcmd(mkdb1)
    runcmd(mkdb2)
    runcmd(blast1)
    runcmd(blast2)


def parse_blast(blast_out):
    my_dt = {}
    blast_num = 0

    with open(blast_out, 'r') as lines:
        for line in lines:
            blast_num += 1
            line = line.strip()
            line_spl = line.split("\t")
            qseqid = line_spl[0]
            sseqid = line_spl[1]
            e_value = float(line_spl[10])

            try:
                my_dt[qseqid]
            except KeyError:
                my_dt[qseqid] = [e_value, sseqid]
            else:
                if my_dt[qseqid][0] > e_value:
                    my_dt[qseqid] = [e_value, sseqid]
#                elif my_dt[qseqid][0] == e_value:
 #                   my_dt[qseqid].append(sseqid)
    return(my_dt, blast_num)

def parse_orthlogs():
    seq1_dt, blast_num1 = parse_blast("seq_1.out")
    seq2_dt, blast_num2 = parse_blast("seq_2.out")
    orthlogs_lst = []
    seq1_lst = []
    seq2_lst = []
    for k,v in seq1_dt.items():
        for i in range(1, len(v)):
            seq1_lst.append(k + "_SeparateMarker_" + v[i])

    for k,v in seq2_dt.items():
        for i in range(1, len(v)):
            seq2_lst.append(v[i] + "_SeparateMarker_" + k)

    for seq in seq1_lst:
        if seq in seq2_lst:
            pairs = seq.replace("_SeparateMarker_", "\t")
            orthlogs_lst.append(pairs)

    r = open("README.txt", "w")
    r.write("the number of initial blast hits is: " + str(int(blast_num1) + int(blast_num2)) + "\n")
    r.write("the number of orthologous gene pairs is: " + str(len(orthlogs_lst)))
    r.close()

    f = open(out_file, "w")
    for lst in orthlogs_lst:
        f.write(lst + "\n")
    f.close()

    subprocess.run("rm seq_[1-2].db*", shell=True, encoding="utf-8")
 #   subprocess.run("rm seq_[1-2].out", shell=True, encoding="utf-8")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'find_orthologs')
    parser.add_argument('--i1', dest = 'i1', help = '<input file 1>')
    parser.add_argument('--i2', dest = 'i2', help = '<input file 2>')
    parser.add_argument('--o', dest = 'o', help = '<Output file name>')
    parser.add_argument('--t', dest = 't', help = '<Sequence type -n/p>')
    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    out_file = args.o
    perform_blast(args)
    parse_orthlogs()

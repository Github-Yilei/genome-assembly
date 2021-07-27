import sys

file = sys.argv[1]
SNPs = int()
Insertions = int()
Deletions = int()

with open(file, 'r') as vcf:
        for line in vcf:
                if line.startswith('#'):
                        continue
                else:
                        line_spl = line.split()
                        ref = len(line_spl[3])
                        alt = len(line_spl[4])
                        if ref == alt:
                                SNPs += 1
                        elif ref < alt:
                                Insertions += 1
                        elif ref > alt:
                                Deletions += 1

print(str(SNPs) + "\t" + str(Insertions) + "\t" + str(Deletions))

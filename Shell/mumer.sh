
Mummer=~/miniconda3/bin
${Mummer}/nucmer --mum -l 1000 -c 200 -g 200  -p hap2 monoploid.fa hap2.fa

# 长度在1000，相似度大于90
# -1: 1对1联配，允许重排，是-r和-q的交集

${Mummer}/delta-filter -1 -i 90 -l 1000 hap2.delta > hap2_filtered.delta
${Mummer}/show-coords -T -q -H -r hap2_filtered.delta > hap2_filtered.coords

# dot plot
#mumer_dot.R

# snp vcf
${Mummer}/show-snps -C -T -l -r ref_qry_filtered.delta > ref_qry_filtered.snps
python3 mummer2vcf.py --input-header -s ref_qry_filtered.snp -g reference.fa --output-header > qry_snp.vcf

# python3 SNP_calculater.py 

for((K=1;K<=9;K++));
do echo "/share/home/stu_wuyilei/miniconda3/bin/nucmer --mum -l 1000 -c 200 -g 200  -p Chr${K} Hap1Chr/chr${K}_RagTag.fasta Hap2Chr/chr${K}_RagTag.fasta && show-snps -Clr Chr${K}.delta > Chr${K}.snps" >> cmd.list; 
done

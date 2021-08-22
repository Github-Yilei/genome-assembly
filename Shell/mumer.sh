
Mummer=~/miniconda3/bin/nucmer

${Mummer}/nucmer -p ${sampleID} -t ${RequiredCPU} reference.fa query.fa
${Mummer}/delta-filter -1 -q -r  ref_qry.delta > ref_qry_filtered.delta
${Mummer}/show-coords -q -r -T -d -H -c ref_qry_filtered.delta > ref_qry_filtered.coords

# dot plot
#mumer_dot.R

# snp vcf
${Mummer}/show-snps -C -T -l -r ref_qry_filtered.delta > ref_qry_filtered.snps
python3 mummer2vcf.py --input-header -s ref_qry_filtered.snp -g reference.fa --output-header > qry_snp.vcf

# python3 SNP_calculater.py 

for((K=1;K<=9;K++));
do echo "/share/home/stu_wuyilei/miniconda3/bin/nucmer --mum -l 1000 -c 200 -g 200  -p Chr${K} Hap1Chr/chr${K}_RagTag.fasta Hap2Chr/chr${K}_RagTag.fasta && show-snps -Clr Chr${K}.delta > Chr${K}.snps" >> cmd.list; 
done

Mummer=/share/home/stu_wuyilei/biosoft/ppsPCP_file/mummer-4.0.0beta2
RequiredCPU=5


${Mummer}/nucmer -p ${sampleID} -t ${RequiredCPU} reference.fa query.fa
${Mummer}/delta-filter -1 -q -r  ref_qry.delta > ref_qry_filtered.delta
${Mummer}/show-coords -q -r -T ref_qry_filtered.delta > ref_qry_filtered.coords

# dot plot
mumer_dot.R

# snp vcf
${Mummer}/show-snps -C -T -l -r ref_qry_filtered.delta > ref_qry_filtered.snps
python3 mummer2vcf.py --input-header -s ref_qry_filtered.snp -g reference.fa --output-header > qry_snp.vcf

python3 SNP_calculater.py 

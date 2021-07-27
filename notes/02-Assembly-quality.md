## mapping rate

The Illumina short reads were mapped back to the assembled contigs using Burrows–Wheeler Aligner software. The highly mapping rate indicating high completeness and accuracy of the final assembly.


## QUAST

```
~/biosoft/quast-5.0.2/quast.py -t 5 -o quast_out ${prefix}.contigs.fasta
```

## busco

To assess the completeness of the assembly, we performed benchmarking universal single-copy ortholog (BUSCO) analysis by searching against the conserved single-copy
genes in plants and identified the number of complete BUSCOs.

```
nohup ~/miniconda3/envs/busco/bin/busco -m genome -l embryophyta_odb10 -o busco_out -c 5 --offline  -i purged.fa &

mkdir summaries
ln -s short_summary_prefix.txt summaries
python ~/miniconda3/envs/busco/bin/generate_plot.py -wd summaries --no_r
# manually editting R scripts
```
## LAI 

LTR assembly index assess assembly continuity, please keep in mind **Not all** species has abundant and avtive LTR retrotransposons to perform LAI. 

```
ltr_finder groups.asm.fasta >groups.asm.finder.scn
				
LTR_retriever/LTR_retriever -threads 4 \
	-genome groups.asm.fasta \
	-infinder groups.asm.finder.scn

LTR_retriever/LAI -t 10 \
	-genome groups.asm.fasta \
	-intact groups.asm.fasta.pass.list \
	-all groups.asm.fasta.out

# firts line is whole genome stats
less groups.asm.fasta.mod.out.LAI
```

## mummer SNPs
To investigate sequence divergence and evolutionary brelationships, we stringently aligned genome sequences with no gaps or indels allowed within an alignment block.

```
Mummer=/share/home/stu_wuyilei/biosoft/ppsPCP_file/mummer-4.0.0beta2
RequiredCPU=5


${Mummer}/nucmer -p ${sampleID} -t ${RequiredCPU} reference.fa query.fa
${Mummer}/delta-filter -1 -q -r  ref_qry.delta > ref_qry_filtered.delta
${Mummer}/show-coords ref_qry_filtered.delta > ref_qry_filtered.coords

# dot plot
mumer_dot.R

# snp vcf
${Mummer}/show-snps -C -T -l -r ref_qry_filtered.delta > ref_qry_filtered.snps
python3 mummer2vcf.py --input-header -s ref_qry_filtered.snp -g reference.fa --output-header > qry_snp.vcf

python3 SNP_calculater.py 
```



## switch errors in the phased genome assembly

A switch error indicates that a single base that is supposed to be present in one haplotype is incorrectly anchored onto another. This kind of assembly error is likely prevalent in the haplotype-resolved genome assembly.

### SNP phasing for PacBio HiFi reads

```
~/miniconda3/bin/minimap2 -t 5 --secondary=no -ax map-hifi hap1_chrom.fa unmapped_combined_CCS.fq.gz -o hap1.sam 

```
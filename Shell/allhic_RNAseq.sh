#!/bin/bash
usage()
{
	echo "    Usage: `basename $0` -g ref.fa -refc ref_cds.fa -refg ref_gene.gff3  -R1 RNA_1.fq -R2 RNA_2.fq -h1 HiC_1.fa -h2 HiC_2.fa"
	echo "          -g: target genome"
  echo "          -refc: CDS of close related species"
	echo "          -refg: gff3 of close related species"
	echo "          -R1: RNA_1.fq"
  echo "          -R2: RNA_2.fq"
  echo "          -h2: HiC reads "
  echo "          -h2: HiC reads"
	echo "          -k: group_count"
	echo "          -e: enzyme_sites (HindIII: AAGCTT; MboI: GATC), default: HindIII"
	echo "          -t: threads, default: 10"
	echo "          -b: bin_size for hic heatmap, can be divided with comma, default: 500k"
	exit 0
}

### get options
while getopts ':r:1:2:k:e:t:b:' OPT; do
	case $OPT in
		g)
			genome="$OPTARG";;
		refc)
			refc="$OPTARG";;
		refg)
			refg="$OPTARG";;
		R1)
			R1="$OPTARG";;
		R2)
			R2="$OPTARG";;
    h1)
			h1="$OPTARG";;
    h2)
			h2="$OPTARG";;
		t)
			threads="$OPTARG";;
		b)
			bin_size="$OPTARG";;
		?)
			usage;;
	esac
done


############## link required files##################
mkdir RNAseq && cd RNAseq
ln -s ${genome} ./draft.asm.fa
ln -s ${refc} ./ref.fa
ln -s ${refg} ./ref_gene.gff
ln -s ${R1} ./RNAseq_R1.fq
ln -s ${R2} ./RNAseq_R2.fq

##############Allele.ctg.table ##################
ref.fa
ref_gene.gff

# build index
~/miniconda3/bin/STAR\
  --runThreadN 30 \
  --runMode genomeGenerate \
  --genomeDir STAR \
  --genomeFastaFiles draft.asm.fa

# alignment
~/miniconda3/bin/STAR\
    --genomeDir STAR \
    --runThreadN 20 \
    --readFilesIn  RNAseq_R1.fq  RNAseq_R2.fq \
    --outFileNamePrefix sample\
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 10 \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --limitBAMsortRAM 40000000000000

# Spliced transcript
~/miniconda3/bin/stringtie sampleAligned.sortedByCoord.out.bam -p 10 -o qry.gtf

# GTF to GFF3
~/miniconda3/pkgs/cufflinks-2.2.1-liulab/bin/gffread qry.gtf -o qry.gff

# get cdna
~/miniconda3/pkgs/cufflinks-2.2.1-liulab/bin/gffread \
        qry.gff -g draft.asm.fa -w qry.fa

# uniform cds name and gene name for classify.pl
sed -e 's/transcript/gene/' -e 's/ID/Name/' qry.gff > qry_gene.gff

# build index
~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/makeblastdb \
        -in ref.fa -dbtype nucl

# alignment
~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/blastn \
        -query qry.fa -db ref.fa \
        -out qry_vs_ref.blast.out -evalue 0.001 -outfmt 6 -num_threads 4 -num_alignments 1

# filter
~/biosoft/ALLHiC/scripts/blastn_parse.pl \
        -i qry_vs_ref.blast.out -o Eblast.out -q qry.fa -b 1 -c 0.6 -d 0.8

# Classify alleles based on BLAST results
~/biosoft/ALLHiC/scripts/classify.pl \
        -i Eblast.out -p 4 -r ref_gene.gff -g qry_gene.gff


############## link required files##################
mv Allele.ctg.table ../
cd ../

ln -s ${genome} ./draft.asm.fasta
ln -s ${h1} ./hic_r1.fq.gz
ln -s ${h2} ./hic_r2.fq.gz

##############allhic ##################
# construct index
~/biosoft/bwa/bwa index -a bwtsw draft.asm.fasta
~/miniconda3/bin/samtools faidx draft.asm.fasta

# alegment
## -5SP
~/biosoft/bwa/bwa mem -t 20 draft.asm.fasta hic_r1.fq.gz hic_r2.fq.gz -o bwa_mem.sam

# filtering sam
~/biosoft/ALLHiC/scripts/PreprocessSAMs.pl bwa_mem.sam draft.asm.fasta HINDIII
~/miniconda3/bin/samtools view -t -b -h -@ 10 bwa_mem.REduced.paired_only.bam -o sampe.clean.sam
~/miniconda3/bin/samtools view -b -t draft.asm.fasta.fai sampe.clean.sam -o sampe.clean.bam

# Prune
~/biosoft/ALLHiC/bin/ALLHiC_prune -i Allele.ctg.table -b sampe.clean.bam -r draft.asm.fasta

# Partition
~/biosoft/ALLHiC/bin/ALLHiC_partition -b prunning.bam -r draft.asm.fasta -e AAGCTT -k 9


# Rescue
~/biosoft/ALLHiC/bin/ALLHiC_rescue -b sampe.clean.bam -r draft.asm.fasta \
        -c prunning.clusters.txt \
        -i prunning.counts_AAGCTT.txt

### optimize
~/biosoft/ALLHiC/bin/allhic extract sampe.clean.bam draft.asm.fasta --RE AAGCTT

rm cmd.list
for((K=1;K<=9;K++))
do echo "~/biosoft/ALLHiC/bin/allhic optimize prunning.counts_AAGCTT.9g${K}.txt sampe.clean.clm" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU 9

### Build
~/biosoft/ALLHiC/bin/ALLHiC_build draft.asm.fasta

### plot
~/miniconda3/bin/samtools faidx groups.asm.fasta
cut -f1,2 groups.asm.fasta.fai| grep 'AAGCTT' > chrn.list

~/biosoft/ALLHiC/bin/ALLHiC_plot sampe.clean.bam groups.agp chrn.list 100k pdf



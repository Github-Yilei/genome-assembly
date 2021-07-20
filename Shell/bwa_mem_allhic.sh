#!/bin/bash

### index reference genome
~/biosoft/bwa/bwa index draft.asm.fasta 
~/miniconda3/bin/samtools faidx draft.asm.fasta 

### bwa mem
# -5SP
~/biosoft/bwa/bwa mem -t 10 draft.asm.fasta HiC_R1.fastq.gz HiC_R2.fastq.gz -o bwa_mem.sam 

### building atg_tabel
~/miniconda3/envs/allhic_env/bin/gmap_build -D . -d DB draft.asm.fasta 
~/miniconda3/envs/allhic_env/bin/gmap -D . -d DB -t 12 -f 2 -n 2 reference.cds.fasta > gmap.gff3
perl gmap2AlleleTable.pl ref_gene.gff3

### filtering sam
~/biosoft/ALLHiC/scripts/PreprocessSAMs.pl bwa_mem.sam draft.asm.fasta HINDIII
~/miniconda3/bin/samtools view -t -b -@ 10 bwa_mem.REduced.paired_only.bam -o sampe.clean.sam
~/miniconda3/bin/samtools view -b -t draft.asm.fasta.fai sampe.clean.sam -o sampe.clean.bam

### Prune	
~/biosoft/ALLHiC/bin/ALLHiC_prune -i Allele.ctg.table -b sampe.clean.bam -r draft.asm.fasta

### Partition
~/biosoft/ALLHiC/bin/ALLHiC_partition -b prunning.bam -r draft.asm.fasta -e AAGCTT -k 9


### Rescue
~/biosoft/ALLHiC/bin/ALLHiC_rescue -b sampe.clean.bam -r draft.asm.fasta \
        -c prunning.clusters.txt \
        -i prunning.counts_AAGCTT.txt

### optimize
~/biosoft/ALLHiC/bin/allhic extract sampe.clean.bam draft.asm.fasta --RE AAGCTT

rm cmd.list
for((K=1;K<=9;K++))
do echo "~/biosoft/ALLHiC/bin/allhic optimize sampe.clean.counts_AAGCTT.9g${K}.txt sampe.clean.clm" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  9

### Build
~/biosoft/ALLHiC/bin/ALLHiC_build draft.asm.fasta

### plot
~/miniconda3/bin/samtools faidx groups.asm.fasta
cut -f1,2 groups.asm.fasta.fai| grep 'AAGCTT' > chrn.list

~/biosoft/ALLHiC/bin/ALLHiC_plot sampe.clean.bam groups.agp chrn.list 500k pdf

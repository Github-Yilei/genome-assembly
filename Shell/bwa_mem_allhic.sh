#!/bin/bash
~/miniconda3/bin/samtools faidx draft.asm.fasta 

# building atg_tabel
~/miniconda3/envs/allhic_env/bin/gmap_build -D . -d DB draft.asm.fasta 
~/miniconda3/envs/allhic_env/bin/gmap -D . -d DB -t 12 -f 2 -n 2 reference.cds.fasta > gmap.gff3
perl gmap2AlleleTable.pl ref_gene.gff3

# filtering sam
~/biosoft/ALLHiC/scripts/PreprocessSAMs.pl bwa_mem.sam draft.asm.fasta MBOI
~/miniconda3/bin/samtools view -t -b bwa_mem.REduced.paired_only.bam > sampe.clean.sam
~/miniconda3/bin/samtools view -b -t draft.asm.fasta.fai sample.clean.sam > sample.clean.bam

# Partition
~/biosoft/ALLHiC/bin/ALLHiC_partition -b sampe.clean.bam -r draft.asm.fasta -e AAGCTT -k 9 

# Optimize 
~/biosoft/ALLHiC/bin/allhic extract sampe.clean.bam draft.asm.fasta --RE AAGCTT 

rm cmd.list
for i in group*.txt; 
do echo "~/biosoft/ALLHiC/bin/allhic optimize $i sampe.clean.clm" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  9

# Build
~/biosoft/ALLHiC/bin/ALLHiC_build draft.asm.fasta

# plot
~/miniconda3/bin/samtools faidx groups.asm.fasta
cut -f 1,2 groups.asm.fasta.fai  >chrn.list

~/biosoft/ALLHiC/bin/ALLHiC_plot sampe.bwa_mem.REduced.paired_only.bam groups.agp chrn.list 500k pdf
~/biosoft/ALLHiC/bin/ALLHiC_plot sampe.clean.bam groups.agp chrn.list 500k pdf

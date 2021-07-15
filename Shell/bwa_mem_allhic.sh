#!/bin/bash
~/miniconda3/bin/samtools faidx draft.asm.fasta 

# building atg_tabel
#~/miniconda3/envs/allhic_env/bin/gmap_build -D . -d DB draft.asm.fasta 
#~/miniconda3/envs/allhic_env/bin/gmap -D . -d DB -t 12 -f 2 -n 2 reference.cds.fasta > gmap.gff3
#perl gmap2AlleleTable.pl ref_gene.gff3

# filtering sam
~/biosoft/ALLHiC/scripts/PreprocessSAMs.pl bwa_mem.sam draft.asm.fasta HINDIII
~/miniconda3/bin/samtools view -t -b -@ 10 bwa_mem.REduced.paired_only.bam > sampe.clean.sam
~/miniconda3/bin/samtools view -b -t -@ 10 draft.asm.fasta.fai sampe.clean.sam > sampe.clean.bam

# Partition
~/biosoft/ALLHiC/bin/ALLHiC_partition -b sampe.clean.bam -r draft.asm.fasta -e AAGCTT -k 9 

# Optimize 
~/biosoft/ALLHiC/bin/allhic extract sampe.clean.bam draft.asm.fasta --RE AAGCTT


### optimize
rm cmd.list
for((K=1;K<=9;K++))
do echo "~/biosoft/ALLHiC/bin/allhic optimize sampe.clean.counts_AAGCTT.9g${K}.txt sampe.clean.clm" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  9

# Build
~/biosoft/ALLHiC/bin/ALLHiC_build draft.asm.fasta

# plot
~/miniconda3/bin/samtools faidx groups.asm.fasta
cut -f1,2 groups.asm.fasta.fai| grep AAGCTT > chrn.list

~/biosoft/ALLHiC/bin/ALLHiC_plot sampe.clean.bam groups.agp chrn.list 500k pdf

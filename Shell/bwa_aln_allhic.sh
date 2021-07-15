### combine SAM
grep ^@ sampe_part1_aln.sam >head
ls | grep aln.sam | while read id; do grep -v ^@ $id ; done > combined.sam
cat head combined.sam > combined_aln.sam
rm head combined.sam

### Filtering SAM
~/biosoft/ALLHiC/scripts/PreprocessSAMs.pl combined_aln.sam draft.asm.fasta HINDIII
~/biosoft/ALLHiC/scripts/filterBAM_forHiC.pl combined_aln.REduced.paired_only.bam sample.clean.sam
~/miniconda3/bin/samtools view -bt draft.asm.fasta.fai sample.clean.sam > sample.clean.bam

### Prune	
~/biosoft/ALLHiC/bin/ALLHiC_prune -i Allele.ctg.table -b sampe.clean.bam -r draft.asm.fasta

### Partition
~/biosoft/ALLHiC/bin/ALLHiC_partition -b prunning.bam -r draft.asm.fasta -e AAGCTT -k 18

### Rescue
~/biosoft/ALLHiC/bin/ALLHiC_rescue -b sampe.clean.bam -r draft.asm.fasta -c clusters.txt -i counts_RE.txt 

### optimize
~/biosoft/ALLHiC/bin/allhic extract sampe.clean.bam draft.asm.fasta --RE AAGCTT

rm cmd.list
for((K=1;K<=18;K++))
do echo "~/biosoft/ALLHiC/bin/allhic optimize sampe.clean.counts_AAGCTT.18g${K}.txt sampe.clean.clm" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  18

# Build
~/biosoft/ALLHiC/bin/ALLHiC_build draft.asm.fasta

# plot
~/miniconda3/bin/samtools faidx groups.asm.fasta
cut -f1,2 groups.asm.fasta.fai| grep AAGCTT > chrn.list

~/biosoft/ALLHiC/bin/ALLHiC_plot sampe.clean.bam groups.agp chrn.list 500k pdf

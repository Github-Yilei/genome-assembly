### Filtering SAM
#for((P=1;P<=8;P++))
#do
#~/biosoft/ALLHiC/scripts/PreprocessSAMs.pl sampe_part${P}_aln.sam draft.asm.fasta HINDIII
#done

rm cmd.list
for((P=1;P<=8;P++))
do
echo "~/biosoft/ALLHiC/scripts/filterBAM_forHiC.pl sampe_part${P}_aln.REduced.paired_only.bam sampe_part${P}.clean.sam" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  8

rm cmd.list
for((P=1;P<=8;P++))
do
echo "~/miniconda3/bin/samtools view -bt draft.asm.fasta.fai sampe_part${P}.clean.sam >sampe_part${P}.clean.bam" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  8

~/miniconda3/bin/samtools merge sampe_part1.clean.bam sampe_part2.clean.bam sampe_part3.clean.bam sampe_part4.clean.bam sampe_part5.clean.bam \
	sampe_part6.clean.bam sampe_part7.clean.bam sampe_part8.clean.bam -O BAM -@ 10 >sampe.clean.bam

### Prune	
~/biosoft/ALLHiC/bin/ALLHiC_prune -i Allele.ctg.table -b sampe.clean.bam -r draft.asm.fasta

### Partition
~/biosoft/ALLHiC/bin/ALLHiC_partition -b prunning.bam -r draft.asm.fasta -e AAGCTT -k 9

### Rescue
~/biosoft/ALLHiC/bin/ALLHiC_rescue -b sampe.clean.bam -r draft.asm.fasta -c clusters.txt -i counts_RE.txt 

### optimize
~/biosoft/ALLHiC/bin/allhic extract sampe.clean.bam draft.asm.fasta --RE AAGCTT

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

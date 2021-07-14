~/miniconda3/bin/samtools faidx draft.asm.fasta

###
rm cmd.list
for i in {1..8}
do echo "~/biosoft/ALLHiC/scripts/filterBAM_forHiC.pl part${i}_aln.REduced.paired_only.bam sample${i}.clean.sam" >> cmd.list
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  9

###
rm cmd.list
for i in {1..8}
do echo "~/miniconda3/bin/samtools view -bt draft.asm.fasta.fai sample${i}.clean.sam | ~/miniconda3/bin/samtools sort -o sample${i}.clean.sorted.bam"
done
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU  9

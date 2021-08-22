##  A Complete Eukaryotic Annotation Pipeline

1. repeat mask
2. ab initio gene prediction
	- [GeneMark](http://topaz.gatech.edu/Genemark/background.html): GeneMark-ES, GeneMarkS, MetaGeneMark, GeneMark.hmm
	- Augustus: the well-studied close-related species or [retrained ab initio model](http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.Augustus) will be used to perform Gene Prediction.
	- [SNAP](https://github.com/KorfLab/SNAP)
	- [GlimmerHMM](http://ccb.jhu.edu/software/glimmerhmm/): GHMM structure includes introns of each phase, intergenic regions, and four types of exons (initial, internal, final, and single)
	- FGENESH is a commercial software.
3. protein homology detection and intron resolution
	- [GenomeThreader](https://genomethreader.org/)
4. RNA-based 
	- [HISAT2](https://daehwankimlab.github.io/hisat2/) + [StringTie](http://ccb.jhu.edu/software/stringtie/)
	- [PASA](https://github.com/PASApipeline/PASApipeline): alignment of known ESTs, full-length cDNAs, and most recently, Trinity RNA-Seq assemblies to the genome.
5. [EVidenceModeler(EVM)](https://github.com/EVidenceModeler/EVidenceModeler) to compute weighted consensus gene structure annotations based on the above (2, 3, 4)
6. PASA to update the EVM consensus predictions, adding UTR annotations and models for alternatively spliced isoforms(leveraging 4 and 5).
7. limited manual refinement of genome annotations using [Apollo](http://www.gmod.org/wiki/WebApollo_Installation)

## Transposable elements

1. Interspersed repeats and low complexity DNA sequences should be masked (default: replaced by Ns) before genome annotation. By the way, I think EDTA will solve most of our problems about Repetitive sequence.
2.  Researchers may want to use RepeatModeler for de novo annotation of non-LTR elements, and supplement these annotations with [SINEBase](https://sines.eimb.ru/) or Repbase.

```
mkdir EDTA
conda activate EDTA && cd EDTA
ln -s absolut_path_of_draft.genome draft.asm.fasta
# From head to toe
perl ~/miniconda3/envs/EDTA/bin/EDTA.pl -genome draft.asm.fasta -species others --sensitive 1 --anno 1 --evaluate 1 --curatedlib SINEBase.fa -threads 20 &

```
### Outputs

1. groups.asm.fasta.mod.EDTA.TElib.fa: A non-redundant TE library
2. groups.asm.fasta.mod.MAKER.masked: This is a genome file with only long TEs (>=1 kb) being masked. You may use this for de novo gene annotations. 
3. groups.asm.fasta.mod.EDTA.final/groups.asm.fasta.mod.tbl: statistics of repeat region of assembly.
4. groups.asm.fasta.mod.EDTA.final/groups.asm.fasta.mod.RepeatModeler.raw.fa.masked: all the annotated repeats have been masked by repeatmodeler and RepeatMasker.

### Debug 

1. ERROR: [Can not recognize this MSU position in the list!](https://github.com/oushujun/EDTA/issues/123)
2. Error: gt is not found in the genometools path !
```
chmod u+x pathto/miniconda3/bin/genometools
```

### AUGUSTUS

The well-studied close-related species will be used to perform Gene Prediction if we get a  "new" genome that AUGUSTUS has not been trained before. [or you can retrain on their own ab initio AUGUSTUS for another species.](http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.Augustus)

```
#!/bin/bash
mkdir 01-augustsus && cd 01-augustsus

conda activate augustus
# split genome
~/miniconda3/bin/seqkit split -i draft.fa.masked 

# parallel working
ls genome.fa.masked.split | while read name; 
do 
echo "~/miniconda3/envs/augustus/bin/augustus --species=arabidopsis --gff3=on --UTR=on draft.fa.masked.split/$name > $name.gff " >> cmd.list; 
done

cpu=`ls draft.fa.masked.split | wc -l `
~/miniconda3/pkgs/parafly-r2013_01_21-1/bin/ParaFly -c cmd.list -CPU ${cpu}

# ID setting
 for i in {1..9};
 do
 cat chr$i.masked.gff | perl ~/miniconda3/envs/augustus/bin/join_aug_pred.pl | grep -v '^#' >> temp.joined.gff;
 done

# sort gff
~/miniconda3/bin/bedtools sort -i temp.joined.gff > augustsus.gff3
```

### GMES-ES

```
#!/bin/bash
mkdir 02-Gmes && cd 02-Gmes

# prediction
gmes_petap.pl --ES --sequence draft.fa.masked --cores 50

# transform genemark.gtf to .gff3
perl genemark_gtf2gff3.pl genemark.gtf >tmp.genemark.gff3

# sort gff
~/miniconda3/bin/bedtools sort -i tmp.genemark.gff3 > genemark.gff3
```

### Genomethreader

```
#!/bin/bash
mkdir 03-Genomethreader && cd 03-Genomethreader

# tblastn
~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/makeblastdb -in draft.fa.masked -parse_seqids -dbtype nucl -out genome.db

~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/tblastn -query all.pep.fa -out all_pep.blast -db genome.db -outfmt 6 -evalue 1e-5 -num_threads 20 -qcov_hsp_perc 50.0 -num_alignments 5

# extract pprotein id
awk '{print $1}' all_pep.blast | sort | uniq >unique_pep.list

# extact uniq protein 
 ~/miniconda3/bin/seqkit grep -f unique_pep.list all.pep.fa -o seqkit_unique_pep.fa

# homology prediction
~/miniconda3/bin/gth -genomic draft.fa.masked -protein seqkit_unique_pep.fa -intermediate -gff3out -o gth_homology.gff 
```

### HISAT2 + StringTie

```
#!/bin/bash

"""
@ $1 the absolut path of masked genome
@ $2 the absolut path of RNA forward reads
@ $3 the absolut path of RNA reverse reads
"""

mkdir 04-RNAbased && cd 04-RNAbased
ln -s $1 draft.fa.masked
ln -s $2 RNAseq_1.fastq.gz
ln -s $3 RNAseq_2.fastq.gz

#build index
mkdir index
~/miniconda3/pkgs/hisat2-2.2.1-he1b5a44_2/bin/hisat2-build draft.fa.masked index/draft.fa.masked

# alligment
~/miniconda3/pkgs/hisat2-2.2.1-he1b5a44_2/bin/hisat2 --dta -p 20 -x index/draft.fa.masked -1 RNAseq_1.fastq.gz -2 RNAseq_2.fastq.gz -S RNA.sam
~/miniconda3/bin/samtools sort -@ 10 RNAseq.sam -O BAM -o RNAseq.bam 

###################################
# repeat it for every RNA pair reads
#hisat2 --dta -p 20 -x index/draft.fa.masked -1 RNAseq_2_R1.fq.gz -2 RNAseq_2_R1.fq.gz

# merge the reslut
#~/miniconda3/bin/samtools merge -@ 10 RNAseq.bam RNAseq1.bam RNAseq2.bam RNAseq3.bam
###################################

#  predicts transcripts
~/miniconda3/bin/stringtie -p 10 -o RNAseq.gtf RNAseq.bam

# prepar the data
conda activate transdecoder_env

#  get fasta
~/miniconda3/envs/transdecoder_env/bin/gtf_genome_to_cdna_fasta.pl RNAseq.gtf draft.fa.masked > transcripts.fasta

# transform to gff3
~/miniconda3/envs/transdecoder_env/bin/gtf_to_alignment_gff3.pl RNAseq.gtf > transcripts.gff3

## extracting long ORFs
~/miniconda3/pkgs/transdecoder-5.5.0-pl526_2/opt/transdecoder/TransDecoder.LongOrfs -t transcripts.fasta

~/miniconda3/pkgs/transdecoder-5.5.0-pl526_2/opt/transdecoder/TransDecoder.Predict -t transcripts.fasta

~/miniconda3/envs/transdecoder_env/bin/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transdecoder.gff3
   
echo all jobs Done at `date`
echo transdecoder.gff3 is the evidence for EvidenceModeler

```

## Basic PASA
```
mkdir 05-PASA && cd 05-PASA

ln -s $1 draft.fa.masked
ln -s $2 RNAseq_1.fastq.gz
ln -s $3 RNAseq_2.fastq.gz

conda activate PASA_env
# Trinity assemble
~/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/Trinity \
--seqType fq --CPU 20 --max_memory 64G \
--left RNAseq1_R1.fq.gz \
--right RNAseq2_R2.fq.gz

# modify "pasa.alignAssembly.Template.txt"
cp ~/biosoft/PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt ./alignAssembly.config
vim alignAssembly.config

#basic PASA 
~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g draft.fa.masked -t Trinity.fasta --ALIGNERS blat

echo database.sqlite.pasa_assemblies.gff3 is the evidence for EvidenceModeler.
```


## Combine the results

```
mkdir 06-EVM && cd 06-EVM

ln -s ../01-augustsus/augustsus.gff3 
ln -s ../02-Gmes/genemark.gff3
ln -s ../03-Genomethreader/genomethreader.gff3
ln -s ../04-RNAbased/transdecoder.gff3
ln -s ../05-PASA/database.sqlite.pasa_assemblies.gff3 PASA.gff3
ln -s ../07-Update/compreh_init_build/compreh_init_build.gff3

# set weights.txt according to /EVidenceModeler1/simple_example
vim weights.txt

# split data
~/miniconda3/envs/EVidenceModeler_env/opt/evidencemodeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
	--genome chrom.fa.masked \
	--gene_predictions augustsus.gff3 \
	--gene_predictions gmes.gff3 \
	--protein_alignments genomethreader.gff3 \
	--transcript_alignments transdecoder.gff3 \
	--transcript_alignments  PASA.gff3 \
	--transcript_alignments  compreh_init_build.gff3 \
	--segmentSize 100000 \
	--overlapSize 10000 \
	--partition_listing partitions_list.out

# build perform code
~/miniconda3/envs/EVidenceModeler_env/opt/evidencemodeler-1.1.1/EvmUtils/write_EVM_commands.pl \
	--genome chrom.fa.masked \
	--weights /share/home/stu_wuyilei/project/Geome_assembel/Papeda/05-GenomeAnnotation/06-EVM/weights.txt \
	--gene_predictions augustsus.gff3 \
	--gene_predictions gmes.gff3 \
	--protein_alignments genomethreader.gff3 \
	--transcript_alignments transdecoder.gff3 \
	--transcript_alignments  PASA.gff3 \
	--transcript_alignments  compreh_init_build.gff3 \
	--output_file_name evm.out \
	--partitions partitions_list.out >commands.list

# running 
~/miniconda3/pkgs/parallel-20170422-pl5.22.0_0/bin/parallel --jobs 20 < commands.list

# combine reslut 
~/miniconda3/envs/EVidenceModeler_env/opt/evidencemodeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl \
	--partitions partitions_list.out \
	--output_file_name evm.out

# transform file to .gff3
~/miniconda3/envs/EVidenceModeler_env/opt/evidencemodeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl  \
	--partitions partitions_list.out \
	--output evm.out  \
	--genome chrom.fa.masked

# set name
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVidenceModeler_combined.gff3
```


## Filter sequence

过滤出长度不足50AA的序列，基于这些序列过滤原来的的注释

```bash
~/miniconda3/envs/cufflinks_env/bin/gffread EVM.all.gff3 -g groups.asm.fasta -y tr_cds.fa

~/miniconda3/bin/bioawk -c fastx 'length($seq)<50 {print $name}' tr_cds.fa | cut -d '=' -f 2 > short_aa_gene_list.txt

grep -v -w -f short_aa_gene_list.txt  EVM.all.gff3 > filter.gff
```

## Statistics

### 1. avarage length

```bash
cat EVM_exon.gff3 | awk '{length += ($5 - $4)} END { print "average = " length/NR}'
```

## Setting genome order

python gffrename.py  EVM_output.gff  prefix  > renamed.gff

最好进行编号改变

```bash
#!/usr/bin/env python3
import re
import sys

if len(sys.argv) < 3:
    sys.exit()

gff = open(sys.argv[1])
prf = sys.argv[2]

count = 0
mRNA  = 0
cds   = 0
exon  = 0

print("##gff-version 3.2.1")
for line in gff:
    if not line.startswith("\n"):
        records = line.split("\t")
        records[1] = "."
    if re.search(r"\tgene\t", line):
        count = count + 10
        mRNA  = 0
         chr_num = records[0]
         # the number of the genes
        gene_id = prf + chr_num + str(count).zfill(7) 
        records[8] = "ID={}".format(gene_id)
    elif re.search(r"\tmRNA\t", line):
        cds   = 0
        exon  = 0
        mRNA  = mRNA + 1
        mRNA_id    = gene_id + "." + str(mRNA)
        records[8] = "ID={};Parent={}".format(mRNA_id, gene_id)
    elif re.search(r"\texon\t", line):
        exon     = exon + 1
        exon_id  = mRNA_id + "_exon_" + str(exon)
        records[8] = "ID={};Parent={}".format(exon_id, mRNA_id)
    elif re.search(r"\tCDS\t", line):
        cds     = cds + 1
        cds_id  = mRNA_id + "_cds_" + str(cds)
        records[8] = "ID={};Parent={}".format(cds_id, mRNA_id)
    else:
        continue

    print("\t".join(records))

gff.close()
```

## Cleanning row data

```bash
# fliter gff3 data - options
~/miniconda3/bin/bioawk -c fastx 'length($seq)>50 {print ">"$name "\t"$comment; print $seq}' tr_cds.fa 


conda activate cufflinks_env
# extract protein sequence(automaticlly choose longest)
~/miniconda3/envs/cufflinks_env/bin/gffread EVM.all.gff -g unmasked.fa -y tr_cds.fa
# 提取CDS序列
gffread EVM.all.gff -g unmasked.fa -x cds.fa
# 获得外显子序列
gffread EVM.all.gff -g unmasked.fa -w exons.fa

# remove dot and duplicate >id id
python remove_dot_Proteins.py kiwi_v1.protein.fa > clean_kiwi_v1.1.protein.fa
```

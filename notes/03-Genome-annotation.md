## Transposable elements

Interspersed repeats and low complexity DNA sequences should be masked (default: replaced by Ns) before genome annotation. By the way, I think EDTA will solve most of our problems about Repetitive sequence.

Data: groups.asm.fasta

```
# From head to toe
perl ~/miniconda3/envs/EDTA/bin/EDTA.pl -genome groups.asm.fasta -species others --sensitive 1 --anno 1 --evaluate 1 -threads 20 &

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

### AUGUSTUS

The well-studied close-related species will be used to perform Gene Prediction if we get a  "new" genome that AUGUSTUS has not been trained before. [or you can retrain on their own ab initio AUGUSTUS for another species.](http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.Augustus)

```
mkdir 01-augustsus && cd 01-augustsus

conda activate augustus_env

# split genome >>> genome.fa.masked.split
~/miniconda3/bin/seqkit split -i genome.fa.masked 

# parallel working
find genome.fa.masked.split/ -type f -name "*.fa" | ~/miniconda3/pkgs/parallel-20170422-pl5.22.0_0/bin/parallel -j 30 augustus --species=arabidopsis --gff3=on >> temp.gff 

# ID setting
join_aug_pred.pl < temp.gff  | grep -v '^#' > temp.joined.gff

~/miniconda3/bin/bedtools sort -i temp.joined.gff > augustsus.gff3
```

### GMES-ES

```
# prediction
gmes_petap.pl --ES --sequence genome.fa.masked --cores 50

# transform genemark.gtf to .gff3
perl genemark_gtf2gff3.pl genemark.gtf >genemark.gff3
```

### Genomethreader
```
genomethreader
gunzip *.pep.all.fa.gz
cat *.pep.all.fa > all.pep.fa

# blast >>> tblastn
~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/makeblastdb -in groups.asm.fasta.masked -parse_seqids -dbtype nucl -out groups.masked &

~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/tblastn -query all.pep.fa -out all_pep.blast -db masked -outfmt 6 -evalue 1e-5 -num_threads 20 -qcov_hsp_perc 50.0 -num_alignments 5

# extract pprotein id
awk '{print $1}' all_pep.blast > all_pep.list
sort all_pep.list | uniq >unique_pep.list

# protein 单行显示
~/miniconda3/bin/seqkit seq all.pep.fa -w 0 > all.fa

# 提取匹配蛋白
cat unique_pep.list | while read line;  do grep $line -A 1 all.fa ; done >unique_pep.fa

# homology prediction
~/miniconda3/bin/gth -genomic groups.asm.fasta.masked -protein unique_pep.fa -intermediate -gff3out > gth_homology.gff
```
*Prepare tabel*

1. 去掉#开头，去掉含"*prime_cis_splice_site"的行，去掉gene后面的Target
2. 同基因组下外显子和内含子需要编号，便于识别。

```
# check types
awk '!/^[#]/ {print $3}' gth_homology.gff | sort | uniq

python genomeThreader_to_evm_gff3.py --input_file gth_homo.gff --output gth_cleaned_homology.gff
```
### HISAT2 + StringTie

```
#build index
mkdir index
~/miniconda3/pkgs/hisat2-2.2.1-he1b5a44_2/bin/hisat2-build \
draft.fa.masked \
index/draft.fa.masked

# alligment
~/miniconda3/pkgs/hisat2-2.2.1-he1b5a44_2/bin/hisat2 --dta -p 20 -x index/groups.asm.fasta.masked -1 ../RNAseq${i}_1.fastq.gz -2 ../RNAseq{$i}_2.fastq.gz | ~/miniconda3/bin/samtools  sort -@ 10 > RNAseq{$i}.bam &

# repeat it for every RNA pair reads
hisat2 --dta -p 20 -x index/draft.fa.masked -1 RNAseq_2_R1.fq.gz -2 RNAseq_2_R1.fq.gz

# merge the reslut
 ~/miniconda3/bin/samtools merge -@ 10 rna_merged.bam \
 RNAseq1.bam RNAseq2.bam RNAseq3.bam

#  predicts transcripts
~/miniconda3/bin/stringtie -p 10 -o merged.gtf merged.bam

# prepar the data
conda active transdecoder_env

#  get fasta
~/miniconda3/envs/transdecoder_env/bin/gtf_genome_to_cdna_fasta.pl merged.gtf ../groups.asm.fasta.masked >transcripts.fasta

# transform to gff3
~/miniconda3/envs/transdecoder_env/bin/gtf_to_alignment_gff3.pl RNA.gtf > transcripts.gff3

## extracting long ORFs
~/miniconda3/pkgs/transdecoder-5.5.0-pl526_2/opt/transdecoder/TransDecoder.LongOrfs \
-t transcripts.fasta

~/miniconda3/pkgs/transdecoder-5.5.0-pl526_2/opt/transdecoder/TransDecoder.Predict \
-t transcripts.fasta

~/miniconda3/envs/transdecoder_env/bin/cdna_alignment_orf_to_genome_orf.pl \
     transcripts.fasta.transdecoder.gff3 \
     transcripts.gff3 \
     transcripts.fasta > transcripts.fasta.transdecoder.genome.gff3
   
# EvidenceModeler need transcripts.fasta.transdecoder.genome.gff3
```

## Combine the results

```
mkdir EVM && cd EVM

# set weights.txt according to /EVidenceModeler1/simple_example
vim weights.txt
## 第一列为来源类型：ABINITIO_PREDICTION, PROTEIN, TRANSCRIPT
## 第二列对应着gff3文件第二列: augustsus, braker
## 第三列为权重 我觉得根据基因组引导组装的ORF的可信度高于组装后比对，所以得分和PASA差不多一样高。从头预测权重一般都是1，但是BRAKER可信度稍微高一点，可以在2~5之间
/share/home/stu_wuyilei/canu_1.6/9_anotation/05-evi/weights.txt

# split data
~/miniconda3/envs/EVidenceModeler_env/opt/evidencemodeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
--genome ../../groups.asm.fasta \
--gene_predictions ../augustsus.gff3 \
--gene_predictions ../genemark.gff3 \
--protein_alignments ../gth_clean_homology.gff \
--transcript_alignments ../transcripts.fasta.transdecoder.genome.gff3 \
--pasaTerminalExons  ../database.sqlite.pasa_assemblies.gff3 \
--segmentSize 100000 \ # --segmentsSize = 基因平均长度加上2个标准差(or 10k)
--overlapSize 10000 \
--partition_listing partitions_list.out
# 可以只运行一次，后面调整weight

# build perform code 
--terminalExons ../database.sqlite.pasa_assemblies.gff3 \

~/miniconda3/envs/EVidenceModeler_env/opt/evidencemodeler-1.1.1/EvmUtils/write_EVM_commands.pl \
--genome ../../groups.asm.fasta \
--weights ~/canu_1.6/9_anotation/06-evi/weights.txt \
--gene_predictions ../augustsus.gff3 \
--gene_predictions ../genemark.gff3 \
--protein_alignments ../gth_clean_homology.gff \
--transcript_alignments ../transcripts.fasta.transdecoder.genome.gff3 \
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
--genome ../../groups.asm.fasta 

# set name
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
```

## PASA

If you have a highly fragmented draft genome, then you are likely better off performing a genome-free de novo transcriptome assembly.Therefore, we can only conduct a basic PASA.

1. Trinity de novo RNA-Seq assemblies

2. [Trinity genome-guided RNA-Seq assemblies](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Genome-Guided-Trinity-Transcriptome-Assembly) ：HISAT2 - RNA-Seq bwa sorted  read alignment bam >>  

   ```bash
    Trinity --genome_guided_bam rnaseq.bam \
            --genome_guided_max_intron 10000 \
            --max_memory 10G --CPU 10 
   ```

3. [StringTie](https://ccb.jhu.edu/software/stringtie/) transcript structures (ex. RNA.gtf)

**Notice：**

a. 2 and 3 were prepared at the step of HISAT2 + StringTie 

b. Trinity to RNA-Seq samples derived from microbial eukaryotes, the '--jaccard_clip' parameter is essential.

c.  database.sqlite should be removed, otherwise will enconter ERROR.

d. absolute path is essential for sqlite database(ie, ../9_anotation/database.sqlite). otherwise, the MySQL database will be used in default

```bash
# dump the fastq file  
 ~/miniconda3/bin/fastq-dump --gzip --split-3 --defline-seq '@$sn[_$rn]/$ri' SRR_id
 
# perform quality control
~/miniconda3/bin/fastp -i -o RNAseq_R1.fq.gz -I -O RNAseq_R2.fq.gz

# Trinity assemble
~/miniconda3/pkgs/trinity-2.8.5-h8b12597_5/opt/trinity-2.8.5/Trinity

Trinity --seqType fq --CPU 50 --max_memory 64G \
--left RNAseq1_R1.fq.gz,RNAseq2_R1.fq.gz,RNAseq2_R1.fq.gz \
--right RNAseq2_R2.fq.gz,RNAseq2_R2.fq.gz,RNAseq3_R2.fq.gz &

# modify "pasa.alignAssembly.Template.txt"
cp ~/biosoft/PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config
vim alignAssembly.config
## DATABASE=database.sqlite (absolute path)
## validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
## validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=80

############ basic PASA  ########################
~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g chi_unmasked.fa -t Trinity.fasta --ALIGNERS blat
# database.sqlite.pasa_assemblies.gff3 is what we need

# other steps
########## Build a Transcriptome Database ###########
cat Trinity.fasta Trinity.GG.fasta > transcripts.fasta

# extract accession id
~/biosoft/PASApipeline/misc_utilities/accession_extractor.pl < Trinity.fasta > tdn.accs

# run PASA
~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl \
-c alignAssembly.config -t transcripts.fasta \
-C -R -g path/draft.fa.masked \
--ALIGNERS blat \
--CPU 2 --TDN tdn.accs --cufflinks_gtf cufflinks.gtf



# cleaning the transcripts[Optional]
~/biosoft/PASApipeline/bin/seqclean  transcripts.fasta

#  PASA alignment
 ~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s 13 -R -g ../groups.asm.fasta -t Trinity.fasta --ALIGNERS blat --CPU 20

~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ../groups.asm.fasta  -t Trinity.fasta --ALIGNERS blat --CPU 20 

# checking the annotations.gff
~/biosoft/PASApipeline/misc_utilities/pasa_gff3_validator.pl orig_annotations.gff3

# Loadding annotations.gff
~/biosoft/PASApipeline/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config -g ../groups.asm.fasta -P EVM.all.gff3

# updatting annotations.gff
~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl -c annotCompare.config -A -g ../groups.asm.fasta  -t Trinity.fasta(all_transcripts.fasta.clean) --CPU 1

```

[SQLite database is locked](https://github.com/PASApipeline/PASApipeline/issues/74): --CPU 1

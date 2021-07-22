## Overviews

This genome assembly strategy employed both Illumina shortread and PacBio long-read sequencing data, with scaffolding informed by Hi–C chromosomal contact information. BioNano reads are optional at filling gap.

## Tools 

- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [fastp](https://github.com/OpenGene/fastp)
- [multiqc](https://multiqc.info/)
- [trimmomatic](https://github.com/timflutre/trimmomatic)
- [pbh5tools](https://github.com/PacificBiosciences/pbh5tools)
- [R-pbh5](https://github.com/PacificBiosciences/R-pbh5) 
- [stsPlots](https://github.com/PacificBiosciences/stsPlotsPlot)
- [jellyfish](https://github.com/gmarcais/Jellyfish)
- [GCE](https://github.com/fanagislab/GCE)
- [canu](https://github.com/marbl/canu)
- [QUAST](http://quast.sourceforge.net/quast.html)
- [HERA](https://github.com/liangclab/HERA), [PBS version](https://github.com/Github-Yilei/HERA)
- [purge_dups](https://github.com/dfguan/purge_dups)
- [pilon](https://github.com/broadinstitute/pilon)
- [Burrows–WheelerAligner](https://github.com/lh3/bwa)
- [samtools](https://github.com/samtools)
- [STAR](https://github.com/alexdobin/STAR)
- [AllHiC](https://github.com/tangerzhang/ALLHiC)
- [fastool](https://github.com/fstrozzi/Fastool)

## Quality control of sequence

Quality control is a key step in high-throughput sequencing experiments, it is should be the first step unless otherwise noted. Not only will this give you an idea of the quality of your data but it will also clean up and reduce the size of your data, making downstream analysis much easier! The main steps in this process are:

- Removing low quality bases
- Removing low complexity reads
- Remove artifacts (barcodes, adapters, chimeras)

The chloroplast and mitochondrion genome sequence MAY be removed from Pacbio reads according to several close-related species. 

**Data:**

- reads of illumina whole genome sequencing

```
cd ${project}/Genome_assembly/workflow
# 1. Perform fastp quality control for 10 samples
sh perform_fastpQC.sh

# or 2. Perform trimmomatic Trimming 
sh trimmomatic_QC.sh

# Summarising the output
cd ${project}/Genome_assembly/${species}/Source_data/Illumina

python3 parse_fastp_json.py --input_path FastpDir --output ./result.txt
~/miniconda3/bin/multiqc RawFastQC -q -o ./RawMultiqc
~/miniconda3/bin/multiqc CleanFastQC -q -o ./CleanMultiqc
```

**Data:**
- reads of Pacbio long reads
- chloroplast and mitochondrion genome from NCBI.

```
  ~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/makeblastdb -in combined.fa -title combined_DB -dbtype nucl -out conbined_database
  ~/miniconda3/pkgs/blast-2.10.1-pl526he19e7b1_1/bin/blastn -db conbined_database -query pacbio.fa -num_threads 20 -evalue 1e-5 -outfmt 6 -out pacbio_blastn.txt
```

## Genome survey

Using Illumina short reads, k-mer distribution was estimated using jellyfish. The overall characteristics of the genome such as genome size, repeat contents, and
heterozygosity rate were estimated using GCE software. 

**Notice:**

- In an actual genome or finished sequence, a k-mer and its reverse complement are not equivalent, hence using the -C switch does not make sense. 

### jellyfish

**Data:**

- reads of illumina whole genome sequencing

```shell
cd ${ProjectDir}/00-GenomeSurvey
seq_name=sequnence_name

# For zip file
~/miniconda3/pkgs/kmer-jellyfish-2.3.0-hc9558a2_1/bin/jellyfish count -m 17 -o jellyfish_17k -s 100M -t 3 -C <(zcat ${seq_name}_1.fq.gz) <(zcat ${seq_name}_2.fq.gz)

# compute the histogram of the k-mer occurrences
~/miniconda3/pkgs/kmer-jellyfish-2.3.0-hc9558a2_1/bin/jellyfish histo -t 3 jellyfish_17k > jellyfish_17k.histo

#  summary states
~/miniconda3/pkgs/kmer-jellyfish-2.3.0-hc9558a2_1/bin/jellyfish stats jellyfish_17k

# checking  
awk '{b = $1*$2; print b"\t"$1}' jellyfish_17k.hist | sort -n -k 1

# plot and statistics
genome_survey.R

```

### GCE

**Data:**

- reads of illumina whole genome sequencing

**notice:**

For the total kmer number for gce option "-g", and the depth frequency file for gce option "-f":
	less prefix.kmer.freq.stat | grep "#Kmer indivdual number" 
	less prefix.kmer.freq.stat | perl -ne 'next if(/^#/ || /^\s/); print; ' | awk '{print $1"\t"$2}' > prefix.kmer.freq.stat.2colum 

```
~/biosoft/gce-1.0.0/kmerfreq -k 17 -l sample_list -t 10 -p prefix -o 0 ./GCE.lib &> kmer_freq.log

~/biosoft/gce-1.0.0/gce -g `value of code 1` -f prefix.kmer.freq.stat.2colum  >gce.table 2>gce.log

# heterzygous mode
~/biosoft/gce-1.0.0/gce -g `value of code 1` -f prefix.kmer.freq.stat.2colum -c 75 -H 1 >gce2.table 2>gce2.log
```

### genomescope

http://qb.cshl.edu/genomescope/

## Removing mitochondrion and chloroplast genome
mitochondrion and chloroplast genomes were downloaded from the NCBI database and these sequences were used to find read sequences which are similar to the PacBio reads by using GMAP or minimap2 aligner at default setting.

```
~/miniconda3/bin/minimap2 -ax map-hifi MitochondrionChloroplast.fa PacBio_ccs.fastq > minimap2.sam

# secondary=no
~/miniconda3/bin/minimap2 -ax map-hifi --secondary=no MitochondrionChloroplast.fa PacBio_ccs.fastq > minimap2_no_secondary.sam
~/miniconda3/bin/samtools fastq -f 4 minimap2.sam -@ 30 -c 6 > unmapped_minimap2_no_secondary.fq.gz

~/miniconda3/bin/samtools fastq -F 4 minimap2.sam -@ 30 -c 6 > mapped_minimap2_no_secondary.fq.gz

# checking reslut
## key informations will be printed on screen or re-got by
## checking secondary mapped reads
~/miniconda3/bin/samtools flagstat minimap2.sam

## checing mapped and unmapped reads
cat acBio_ccs.fastq | echo $((`wc -l`/4))
cat unmapped_minimap2_no_secondary.fq.gz | echo $((`wc -l`/4))
cat mapped_minimap2_no_secondary.fq.gz | echo $((`wc -l`/4))

## cheching with blastn
~/miniconda3/bin/seqkit fq2fa mapped.fq.gz -o mapped.fa -j 10
nohup /share/home/stu_wuyilei/biosoft/ppsPCP_file/ncbi-blast-2.11.0+/bin/blastn -db combined_database -query mapped.fa -num_threads 10 -evalue 1e-6 -outfmt '6 qseqid sseqid pident nident qlen slen evalue bitscore' -out mapped_blastn.txt > mapped_blastn_log.file  2>&1 &

~/miniconda3/bin/seqkit fq2fa unmapped.fq.gz -o unmapped.fa -j 10
nohup /share/home/stu_wuyilei/biosoft/ppsPCP_file/ncbi-blast-2.11.0+/bin/blastn -db combined_database -query unmapped.fa -num_threads 10 -evalue 1e-6 -outfmt '6 qseqid sseqid pident nident qlen slen evalue bitscore' -out unmapped_blastn.txt >unmapped_blastn_log.file  2>&1 &

# making plot data frame
python blast2plot.py --blast mapped_blastn.txt --out parsed_mapped_blastn.txt
awk '{ print $4}' parsed_mapped_blastn.txt | sort -n
```

## assembley

### HiFiasm

Hifiasm is a fast haplotype-resolved de novo assembler for PacBio HiFi reads. Hifiasm does not perform scaffolding for now. You need to run a standalone scaffolder such as SALSA, AllHic or 3D-DNA to scaffold phased haplotigs.

**Notice**:
- When parental short reads are available, hifiasm can also generate a pair of haplotype-resolved assemblies with trio binning. 
- Hifiasm can generate a pair of haplotype-resolved assemblies with paired-end Hi-C reads or Strand-seq.
- [purge_dup](https://github.com/chhylp123/hifiasm/issues/55#issuecomment-736976724)

```
/share/home/stu_wuyilei/miniconda3/bin/hifiasm -o hifiasm -t 40 PacBio_ccs.fq.gz
# get primary contigs in FASTA
awk '/^S/{print ">"$2;print $3}' hifiasm_ctg.gfa >hifiasm_ctg.fa

blast2plot.R
```

与其他基于图形的汇编程序不同，HiFiasm致力于保持所有单倍型的连续性。 HiCanu只试图保持一个亲本单倍型的连续性，并且经常破坏另一个单倍型的连续性，当分离亲本单倍型时，这些突变点将导致单倍型分解的碎片—HiCanu没有充分利用HiFi Reads Hifiasm针对HiFi特点而开发，在hifi数据的组装表现上较同类软件更为突出，在多个基因组上表现出了更高的准确性和组装的连续性。 

https://github.com/chhylp123/hifiasm/issues/46 

https://github.com/chhylp123/hifiasm#hi-c-integration 

https://github.com/tangerzhang/ALLHiC/issues/86

### canu 

Canu assembles reads from PacBio RS II or Oxford Nanopore MinION instruments into uniquely-assemblable contigs, unitigs.

**Notice**:
- HiCanu has support for PacBio HiFi data by specify -pacbio-hifi.
- HiCanu consensus sequences using PacBio HiFi data are typically well above 99.99% We discourage any post-processing/polishing of these assemblies as mis-mapping within repeats can introduce errors.
- The -p option, to set the file name prefix of intermediate and output files.
- The -d option, to create the assembly-directory and run in that directory.
- The -s option will import a list of parameters from the supplied specification (‘spec’) file that will be applied before any from the command line are used.

**Data:**

- PacBio RS II reads：PacBio.fastq.gz

```
# basic usage
prefix=species_name

~/biosoft/canu-2.1/bin/canu -p ${prefix} -d ${prefix} \
	# gnuplotTested=true \
	  maxThreads=20 \
	  genomeSize=<number>[g|m|k] \
	 -pacbio-raw PacBio.fastq.gz  useGrid=false

# for heterozygous genomes
# built a polyploid pipeline
~/biosoft/canu-2.1/bin/canu batOptions="-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \
	-p ${prefix} -d ${prefix} \
	# gnuplotTested=true \
	maxThreads=20 \
	genomeSize=653m \
	-pacbio-raw PacBio.fastq.gz
```

## Duplication 

it is essential that using purge_dups to remove duplication after canu assembley.

https://github.com/marbl/canu/issues/1814

**Data：**

- PacBio RS II reads：PacBio.fastq.gz
- long contigs: ${prefix}.contigs.fasta

```
########## step a1, a2 can Work in parallel #############
prefix=species_name
# step a1. align pacbio data to contigs
## For a large genome, setting minimap2 -I option
~/biosoft/quast-5.0.2/quast_libs/minimap2/minimap2 -t 20 -x map-pb ${prefix}.contigs.fasta PacBio.fastq.gz | gzip -c - > ${prefix}.paf.gz
~/miniconda3/bin/minimap2 -x map-hifi -t 10 ${prefix}.contigs.fasta PacBio.fastq.gz --secondary=no -o ${prefix}.paf

# pbcstat
~/biosoft/purge_dups/bin/pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
~/biosoft/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log

# step a2.Split an assembly and do a self-self alignment
~/biosoft/purge_dups/bin/split_fa ${prefix}.contigs.fasta > ${prefix}_split
~/biosoft/quast-5.0.2/quast_libs/minimap2/minimap2 -xasm5 -DP ${prefix}_split ${prefix}_split | gzip -c - > ${prefix}.split.self.paf.gz
 python3 ~/biosoft/purge_dups/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png
~/biosoft/purge_dups/bin/calcuts calcuts -l36 -m40 -u240 PB.stat> cutoffs 2>calcults.log

# step 2 Purge haplotigs and overlaps
~/biosoft/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov ${prefix}.split.self.paf.gz > dups.bed 2> purge_dups.log

# step 3 remove haplotypic duplications
~/biosoft/purge_dups/bin/get_seqs -e dups.bed ${prefix}.contigs.fasta 

# purged.fa：clean data
# hap.fa：junk, haplotig, and duplication

# Step 4. Merge hap.fa and $hap_asm and redo the above steps to get a decent haplotig set.
cat *.fasta > total.fasta

## checking 
grep ">" total.fasta
### if '>' had not seperated properlly, add a blank line at the end of fasta file. 
for i in $(ls *.fasta);
do
cat ${i}
echo
done > ./temp/total.fasta

~/biosoft/quast-5.0.2/quast.py -t 20 -o quast_out purged.fa
```


### khaper
# Illumina
# https://github.com/lardo/khaper
~/biosoft/khaper/Bin/Graph.pl pipe -i fq.lst -m 2 -k 15 -s 1,3 -d Kmer_15

### busco

```
nohup ~/miniconda3/envs/busco/bin/busco -m genome -l /share/home/stu_wuyilei/project/Geome_assembel/kiwi/5_evaluate_quality/embryophyta_odb10 -o busco_out -c 5 --offline  -i purged.fa &

mkdir summaries
ln -s short_summary_prefix.txt summaries
python ~/miniconda3/envs/busco/bin/generate_plot.py –wd my_summaries  --no_r
# manually editting R scripts
```
## Fillling gaps

HERA is a local assembly tool using assembled contigs and self-corrected long reads as input to resolves repeats efficiently by constructing a connection graph from an overlap graph and filling gaps.

**Data:**

- PacBio RS II reads：PacBio.fastq.gz
- long contigs: ${prefix}.contigs.fasta

```
cp ~/biosoft/HERAFile/HERA/PBS_pipline.sh ./
sh PBS_pipline.sh

~/biosoft/quast-5.0.2/quast.py -t 20 -o quast_out ${prefix}-Final_Genome_HERA.fasta
```
![](https://user-images.githubusercontent.com/82864917/122715233-29eb9180-d29b-11eb-93a4-3e06bbf44cde.png)


## polish

Correcting contiges based on self-alignment and depth, 1 Mb of genomes corresponds to 1 Gb of memory.

**Data:**

- illumina whole genome sequence
- long contigs: ${prefix}.contigs.fasta


```
# create `ln -s` for contigs and illumina reads
ln -s path_to_${prefix}.contigs.fasta ./draft.asm.fa
ln -s path_to_illumina_reads ./reads_1.fq.gz ./reads_2.fq.gz

# create an index 
~/biosoft/bwa/bwa index draft.asm.fa

# mapping reads to the contig
# -t, -@ is hte number of threads
~/biosoft/bwa/bwa mem -t 20 draft.asm.fa  reads_1.fq.gz reads_2.fq.gz | ~/miniconda3/bin/samtools view -S -b | ~/miniconda3/bin/samtools sort -@ 20 -o mapping_sorted.bam
# creating indexed bam
~/miniconda3/bin/samtools index -@ 10 mapping_sorted.bam


################# if not a PCR-free library, removing PCR duplication
# acturally, it already done at purge_dup
# sambamba markdup -t 10 -r mapping.sorted.bam filterd.bam
## build index
# ~/miniconda3/bin/samtools sort -@ 20 mapping_filtered.bam
# ~/miniconda3/bin/samtools index -@ 10  mapping_filtered.bam
########################################


# call pilon
nohup java -Xmx60G -jar ~/biosoft/pilon-1.23.jar --genome draft.asm.fa --fix all --changes --frags mapping_sorted.bam --threads 25 --outdir ./pilon_out >pilon.log &

# 1GB per megabase
```
You can repeat this for several rounds by using round1.pilon round2.pilon

## Evaluate quality

QUAST is much more informative if at least a close reference genome is provided along with the assemblies.

**Data**:

- long_condraft.asm.fasta

```
# without parameter
~/biosoft/quast-5.0.2/quast.py -t 20 -o quast_out draft.asm.fasta

# with close_species_reference.fasta and pilon data 
~/biosoft/quast-5.0.2/quast.py -t 20 -o quast_pilon \
-r ../A.eriantha_genomic.fna draft.asm.fasta 
pilon_round1.fasta pilon_round2.fasta pilon_round3.fasta pilon_round4.fasta
```


## ALLHic Assmebly

Hi-C is used to assemble genomes by using the 3D proximity relationships between short contigs and scaffolds (draft fragment sequences) to anchor, order and orient them. ALLHiC relies on an allelic contig table (Allele.ctg.table) to remove noisy Hi-C signals. However, if you are working on a **diploid genome**, you can ignore Allele.ctg.table and use the suggested [pipeline](https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-scaffolding-of-a-simple-diploid-genome).


**Data:**
- long_contigs: draft.asm.fasta
- Hi-C reads
- clean RNAseq reads: 1. RNAseq_R1.fq.gz, 2. RNAseq_R2.fq.gz
- CDS from related species: reference.cds.fasta
- gff3 from related species: reference.cds.gff

### Filter low quality(optional)

```
# quality score < 20
# adapter
# shorter than  30bp
~/miniconda3/bin/fastp \
	-i reads_R1.fastq.gz -o reads_R1.fastq.gz \
	-I reads_R2.fastq.gz  -O reads_R2.fastq.gz \
	-q 20 -l 30 -w 10 &
```

### building Allele.ctg.table

There are a couple of ways to generate this table. 

#### based on BLAST

[This a BLAST-based method, which requires a chromosomal level assembly of closely related genome.](https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-identify-allelic-contigs)

```
cd ${Scaffolds}
ln -s ${Contigs}/contigs.fa ./draft.asm.fa
# build index
~/miniconda3/bin/STAR\
  --runThreadN 10 \
  --runMode genomeGenerate \
  --genomeDir STAR \
  --genomeFastaFiles draft.asm.fa
  
# alignment
~/miniconda3/bin/STAR\
    --genomeDir STAR \
    --runThreadN 10 \
    --readFilesIn  RNAseq_R1.fq.gz  RNAseq_R2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix sample\
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 10 \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical

# Spliced transcript
~/miniconda3/bin/stringtie sampleAligned.sortedByCoord.out.bam -p 10 -o qry.gtf

# GTF to GFF3
~/miniconda3/pkgs/cufflinks-2.2.1-liulab/bin/gffread qry.gtf -o qry.gff

# get cdna
~/miniconda3/pkgs/cufflinks-2.2.1-liulab/bin/gffread \
qry.gff -g draft.asm.fa -w qry.fa

# uniform cds name and gene name for classify.pl
sed -e 's/transcript/gene/' -e 's/ID/Name/' qry.gff > qry_gene.gff

# for close related species
ref=ref.gff
grep '[[:blank:]]mRNA[[:blank:]]' $ref | sed -e 's/mRNA/gene/' -e 's/ID/Name/' > ref_gene.gff

sed -e 's/ gene_id "/''/ -e 's/"; gene_version"/''/

ref=ref.gff
grep '[[:blank:]]mRNA[[:blank:]]' $ref | sed -e 's/mRNA/gene/' -e 's/Name/ID/' -e's/ID=transcript:/Name=/' > ref_gene.gff

# now, we get ref.fa, ref_gene.gff, qry.fa, qry_gene.gff
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
```
After running the scripts above, two tables will be generated whose are Allele.gene.table and Allele.ctg.table.

### Perform AllHiC

```
nohup sh bwa_mem_allhic.sh &
```

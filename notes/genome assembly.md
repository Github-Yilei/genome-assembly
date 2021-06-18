## Overviews


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

## Quality control of NGS

Quality control is a key step in high-throughput sequencing experiments, it is should be the first step unless otherwise noted. Not only will this give you an idea of the quality of your data but it will also clean up and reduce the size of your data, making downstream analysis much easier! The main steps in this process are:

- Removing low quality bases
- Removing low complexity reads
- Remove artifacts (barcodes, adapters, chimeras)

**Data:**

- reads of illumina whole genome sequencing

```
# 1. Perform fastp quality control for 10 samples
sh perform_fastpQC.sh 10 

# or 2. Perform trimmomatic Trimming 
sh trimmomatic_QC.sh

# Summarising the output
cd project

python3 parse_fastp_json.py --input_path FastpDir --output ./result.txt
~/miniconda3/bin/multiqc RawFastQC -q -o ./RawMultiqc
~/miniconda3/bin/multiqc CleanFastQC -q -o ./CleanMultiqc
```

## Genome survey

Using Illumina short reads, k-mer distribution was estimated using jellyfish. The overall characteristics of the genome such as genome size, repeat contents, and
heterozygosity rate were estimated using GCE software.

### jellyfish

**Data:**

- reads of illumina whole genome sequencing

```
seq_name=sequnence_name

# For zip file
~/miniconda3/pkgs/kmer-jellyfish-2.3.0-hc9558a2_1/bin/jellyfish count -m 17 -o jellyfish_17k -s 100M -t 3 -c 8 -C <(zcat ${seq_name}_1.fastq.gz) <(zcat ${seq_name}_2.fastq.gz)

# compute the histogram of the k-mer occurrences
~/miniconda3/pkgs/kmer-jellyfish-2.3.0-hc9558a2_1/bin/jellyfish histo -t 3 jellyfish_17k > jellyfish_17k.histo

#  summary states
~/miniconda3/pkgs/kmer-jellyfish-2.3.0-hc9558a2_1/bin/jellyfish stats jellyfish_17k

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

## Working with PACBIO sequencing data

The bas.h5 file and associated bax.h5 files are the main output files produced by the primary analysis pipeline on the PacBio® RS II which contain base-call information. the cmp.h5 is the primary
file is sequence alignment file for SMRT™ sequencing data.

### pbh5tools

pbh5tools is a collection of tools that can manipulate the content or extract data from cmp.h5 or bas.h5:

- bash5tools.py can extract read sequences and quality values for both Raw and circular consensus sequencing (CCS) readtypes and use create `fastq` and `fasta` files.
- cmph5tools.py is a multi-command line tool that can check validity, merge, sort, select, compare, summarize, and stats of cmp.h5 file.

### R-pbh5

R-pbh5 is an R package for parseing HDF5 from the Pacific Biosciences. 

### stsPlots

stsPlots allow the user to assess chip loading, readlength, read score, SNR, and oxygen exclusion to assess potential SMRTcell loading problems. seeing  the associated powerpoint for more information.

## canu assembley

Canu assembles reads from PacBio RS II or Oxford Nanopore MinION instruments into uniquely-assemblable contigs, unitigs.

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
	 -pacbio-raw PacBio.fastq.gz

# for heterozygous genomes
# built a polyploid pipeline
~/biosoft/canu-2.1/bin/canu batOptions="-dg 3 -db 3 -dr 1 -ca 500 -cp 50" \
	-p ${prefix} -d ${prefix} \
	# gnuplotTested=true \
	maxThreads=20 \
	genomeSize=653m \
	-pacbio-raw PacBio.fastq.gz
  
# evaluate assembley
~/biosoft/quast-5.0.2/quast.py -t 20 -o quast_out ${prefix}.contigs.fasta
```

## Duplication 

it is essential that using purge_dups to remove duplication after canu assembley.

**Data：**

- PacBio RS II reads：PacBio.fastq.gz
- long contigs: ${prefix}.contigs.fasta

```
########## step a1, a2 can Work in parallel #############
prefix=species_name
# step a1. align pacbio data to contigs
## For a large genome, setting minimap2 -I option
~/biosoft/quast-5.0.2/quast_libs/minimap2/minimap2 -t 20 -x map-pb ${prefix}.contigs.fasta PacBio.fastq.gz | gzip -c - > ${prefix}.paf.gz

# pbcstat
~/biosoft/purge_dups/bin/pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
~/biosoft/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log

# step a2.Split an assembly and do a self-self alignment
~/biosoft/purge_dups/bin/split_fa ${prefix}.contigs.fasta > ${prefix}_split
~/biosoft/quast-5.0.2/quast_libs/minimap2/minimap2 -xasm5 -DP ${prefix}_split ${prefix}_split | gzip -c - > ${prefix}.split.self.paf.gz

# step 2 Purge haplotigs and overlaps
~/biosoft/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov  ${prefix}.split.self.paf.gz > dups.bed 2> purge_dups.log

# step 3 remove haplotypic duplications
~/biosoft/purge_dups/bin/get_seqs dups.bed ${prefix}.contigs.fasta 

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

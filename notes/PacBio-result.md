## Main output
The bas.h5 file and associated bax.h5 files are the main output files produced by the primary analysis pipeline on the PacBio® RS II which contain base-call information. the cmp.h5 is the primary file is sequence alignment file for SMRT™ sequencing data.

1. pbh5tools

[pbh5tools](https://github.com/PacificBiosciences/pbh5tools) is a collection of tools that can manipulate the content or extract data from cmp.h5 or bas.h5:

bash5tools.py can extract read sequences and quality values for both Raw and circular consensus sequencing (CCS) readtypes and use create fastq and fasta files.
cmph5tools.py is a multi-command line tool that can check validity, merge, sort, select, compare, summarize, and stats of cmp.h5 file.

2. R-pbh5

[R-pbh5](https://github.com/PacificBiosciences/R-pbh5) is an R package for parseing HDF5 from the Pacific Biosciences.

3. stsPlots

[stsPlots](https://github.com/PacificBiosciences/stsPlotsPlot) allow the user to assess chip loading, readlength, read score, SNR, and oxygen exclusion to assess potential SMRTcell loading problems. seeing the associated powerpoint for more information.

## Standard result of subreads.bam

subreads.bam.pbi
subreadset.xml
subreads.bam

bam2fastx, samtools fasta could be used to tranform subreads.bam to fasta or fastaq file.

## Standard result of HiFi

HiFi means highly Accurate Single-Molecule Consensus Reads:

- ccs.len.dis.png
- ccs.stat.xls
- ccs.fastq
- ccs.len.dis.pdf
- ccs.len.dis.xls
- ccs.fasta.len
- ccs.bam

If there is not the CCS files, one can generate HiFi Reads with PacBio [CCS](https://github.com/PacificBiosciences/ccs/blob/develop/docs/how-does-ccs-work.md).

## Uzip and keep raw tar

```
tar -xvf PacBio.tar
```

## Removing mitochondrion and chloroplast genome

mitochondrion and chloroplast genomes were downloaded from the NCBI database and these sequences were used to find read sequences which are similar to the PacBio reads by using GMAP or minimap2 aligner at default setting.

``` 
~/miniconda3/bin/minimap2 -ax map-hifi MitochondrionChloroplast.fa PacBio_ccs.fastq > minimap2.sam
~/miniconda3/bin/samtools fastq -f 4 minimap2.sam -@ 30 -c 6 > unmapped.fq.gz

# key informations will be printed on screen or re-get by
~/miniconda3/bin/samtools flagstat minimap2.sam
cat unmapped.fq.gz | echo $((`wc -l`/4))
```

## Notes

1. HiCanu has support for PacBio HiFi data by specify [-pacbio-hifi](https://canu.readthedocs.io/en/latest/quick-start.html?highlight=hifi#assembling-pacbio-hifi-with-hicanu).

2. HiCanu consensus sequences using PacBio HiFi data are typically well above 99.99% We discourage any post-processing/polishing of these assemblies as mis-mapping within repeats can introduce [errors](https://canu.readthedocs.io/en/latest/quick-start.html?highlight=hifi#consensus-accuracy).

3. HiFiasm
如果没有其他数据，Hifiasm在输出序列时会任意选择每个气泡的一侧输出类似Falcon unzip和HiCanu的主要组装结果（primary contigs）。如果同时有父母本的测序数据，Hifiasm可以通过亲本特有的kmer在图上识别出来自父母本的序列，从而得到两套单倍体基因组。
当然HiFiasm文章中也提到了：

与其他基于图形的汇编程序不同，HiFiasm致力于保持所有单倍型的连续性。
HiCanu只试图保持一个亲本单倍型的连续性，并且经常破坏另一个单倍型的连续性，当分离亲本单倍型时，这些突变点将导致单倍型分解的碎片—HiCanu没有充分利用HiFi Reads
Hifiasm针对HiFi特点而开发，在hifi数据的组装表现上较同类软件更为突出，在多个基因组上表现出了更高的准确性和组装的连续性。
https://github.com/chhylp123/hifiasm/issues/46
https://github.com/chhylp123/hifiasm#hi-c-integration
https://github.com/tangerzhang/ALLHiC/issues/86

## Papers

https://www.pacb.com/wp-content/uploads/Kronenberg-ASHG-2019-High-quality-human-genomes-achieved-through-HiFi-sequence-data-and-FALCON-Unzip-assembly.pdf

https://www.biorxiv.org/content/10.1101/2020.03.14.992248v2

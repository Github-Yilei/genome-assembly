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

## Canu assembly

1. HiCanu has support for PacBio HiFi data by specify [-pacbio-hifi](https://canu.readthedocs.io/en/latest/quick-start.html?highlight=hifi#assembling-pacbio-hifi-with-hicanu).

2. HiCanu consensus sequences using PacBio HiFi data are typically well above 99.99% We discourage any post-processing/polishing of these assemblies as mis-mapping within repeats can introduce [errors](https://canu.readthedocs.io/en/latest/quick-start.html?highlight=hifi#consensus-accuracy).

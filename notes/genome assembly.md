## Overviews


## Working with PACBIO sequencing data

The bas.h5 file and associated bax.h5 files are the main output files produced by the primary analysis pipeline on the PacBio® RS II which contain base-call information. the cmp.h5 is the primary
file is sequence alignment file for SMRT™ sequencing data.

### pbh5tools

[pbh5tools](https://github.com/PacificBiosciences/pbh5tools) is a collection of tools that can manipulate the content or extract data from cmp.h5 or bas.h5:

- bash5tools.py can extract read sequences and quality values for both Raw and circular consensus sequencing (CCS) readtypes and use create `fastq` and `fasta` files.
- cmph5tools.py is a multi-commandline tool that can check validiton, merge, sort, select, compare, summarize and stats cmp.h5

### R-pbh5

[R-pbh5](https://github.com/PacificBiosciences/R-pbh5) is an R package for parseing HDF5 from the Pacific Biosciences. 

### stsPlots

[stsPlots](https://github.com/PacificBiosciences/stsPlotsPlot) allow the user to assess chip loading, readlength, read score, SNR, and oxygen exclusion to assess potential SMRTcell loading problems. seeing  the associated powerpoint for more information.


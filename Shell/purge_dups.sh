#!/bin/bash

usage()
{
	echo "    Usage: `basename $0` draft.asm.fa"
	exit 0
}

prefix=draft
~/miniconda3/bin/minimap2 -x map-hifi -t 10 $1 hifi.fq -o ${prefix}.paf
# pbcstat
~/biosoft/purge_dups/bin/pbcstat *.paf
~/biosoft/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log
python3 ~/biosoft/purge_dups/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png

# step a2.Split an assembly and do a self-self alignment
~/biosoft/purge_dups/bin/split_fa $1 > ${prefix}_split
~/biosoft/quast-5.0.2/quast_libs/minimap2/minimap2 -xasm5 -DP ${prefix}_split ${prefix}_split | gzip -c - > ${prefix}.split.self.paf.gz

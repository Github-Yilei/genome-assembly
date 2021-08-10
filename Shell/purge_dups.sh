#!/bin/bash

usage()
{
        echo "   Usage: `basename $0` -r draft.asm.fa -p draft -h hifi.fq"
        echo "   -r reference genome "
        echo "   -p output PAF to FILE"
        echo "   -h PacBio HiFi reads"
        exit 0
}

### get options
while getopts ':r:p:h:' OPT; do
        case $OPT in
                r)
                        ref="$OPTARG";;
                p)
                        prefix="$OPTARG";;
                h)
                        hifi="$OPTARG";;
                ?)
                        usage;;
        esac
done

# minimap mapping
~/miniconda3/bin/minimap2 -x map-hifi -t 40 ${ref} ${hifi} -o ${prefix}.paf

# pbcstat
~/biosoft/purge_dups/bin/pbcstat *.paf
~/biosoft/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log
python3 ~/biosoft/purge_dups/scripts/hist_plot.py -c cutoffs PB.stat PB.cov.png

# step a2.Split an assembly and do a self-self alignment
~/biosoft/purge_dups/bin/split_fa ${ref} > ${prefix}_split
~/biosoft/quast-5.0.2/quast_libs/minimap2/minimap2 -xasm5 -DP ${prefix}_split ${prefix}_split | gzip -c - > ${prefix}.split.self.paf.gz

# step 2 Purge haplotigs and overlaps
~/biosoft/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov ${prefix}.split.self.paf.gz > dups.bed 2> purge_dups.log

# step 3 remove haplotypic duplications
~/biosoft/purge_dups/bin/get_seqs -e dups.bed ${ref}


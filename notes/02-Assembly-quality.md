First, the Illumina short reads were mapped back to the assembled contigs using Burrows–Wheeler Aligner software. The mapping rate of pairedend reads reached 97.2%, indicating high completeness and accuracy of the final assembly.


## QUAST

```
~/biosoft/quast-5.0.2/quast.py -t 5 -o quast_out ${prefix}.contigs.fasta
```

## busco

```
nohup ~/miniconda3/envs/busco/bin/busco -m genome -l embryophyta_odb10 -o busco_out -c 5 --offline  -i purged.fa &

mkdir summaries
ln -s short_summary_prefix.txt summaries
python ~/miniconda3/envs/busco/bin/generate_plot.py -wd summaries --no_r
# manually editting R scripts
```


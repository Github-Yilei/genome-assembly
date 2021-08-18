
# INFERNAL
```
# Format a CM database
## Pressed and indexed XX CMs and XX HMM filters
~/miniconda3/pkgs/infernal-1.1.3-h516909a_0/bin/cmpress Rfam.cm

## For cmscan 
CMmumber=`cat Rfam.cm | grep 'NAME'| sort | wc -l`
##  Z = (Total ofresidues X 2 X ${CMmumber}/1,000,000)


# NCBI code
~/miniconda3/pkgs/infernal-1.1.3-h516909a_0/bin/cmscan -Z  --nohmmonly --cpu 20  --rfam --cut_ga --fmt 2 --oclan --oskip --clanin Rfam.clanin -o genome.cmscan.out --tblout genome.cmscan.tblout Rfam.cm genome.fa

cmscan -Z 5.874406 --cut_ga --rfam --nohmmonly --tblout genome.tblout --fmt 2 --cpu 8 --clanin Rfam.clanin Rfam.cm genome.fa > genome.cmscan
```

# run PASA
```
~/biosoft/PASApipeline/scripts/Launch_PASA_pipeline.pl \
	-c alignAssembly.config -t transcripts.fasta \
	-C -R -g chrom.fa.masked \
	--ALIGNERS blat \
	--CPU 2 --TDN tdn.accs --cufflinks_gtf cufflinks.gtf
```
## Genomethreader
```
~/biosoft/ppsPCP_file/ncbi-blast-2.11.0+/bin/tblastn -query all.pep.fa -out all_pep.blast -db masked -outfmt 6 -evalue 1e-5 -num_threads 20 -qcov_hsp_perc 50.0 -num_alignments 5
```

BLASTN-RNA

```
~/biosoft/ppsPCP_file/ncbi-blast-2.11.0+/bin/blastn -query ncRNA.fa -db chrom_db -out ncRNA.blast -outfmt 6 -evalue 1e-5
```


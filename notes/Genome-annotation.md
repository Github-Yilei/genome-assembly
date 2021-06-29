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
	- HISAT2 + StringTie
	- PASA: alignment of known ESTs, full-length cDNAs, and most recently, Trinity RNA-Seq assemblies to the genome.
5. EVidenceModeler(EVM) to compute weighted consensus gene structure annotations based on the above (2, 3, 4)
6. PASA to update the EVM consensus predictions, adding UTR annotations and models for alternatively spliced isoforms(leveraging 4 and 5).
7. limited manual refinement of genome annotations using Apollo

## Transposable elements

interspersed repeats and low complexity DNA sequences should be masked (default: replaced by Ns) before genome annotation. By the way, I think EDTA	 will solve most of our problems about Repetitive sequence.

#usr/bin/bash
echo making dir for genome assembly

################################
### Variable information ###
ProjectDir=/share/home/stu_wuyilei/Projects/GenomeAssembly/Species

################################
### Project information ###
Source_data=${ProjectDir}/Source_data

### genome  assembly
GenomeSurvey=${ProjectDir}/00-GenomeSurvey
Contigs=${ProjectDir}/01-Contigs
Scaffolds=${ProjectDir}/02-Scaffolds
PseudoChromosomes=${ProjectDir}/03-PseudoChromosomes

EvaluateQuality=${ProjectDir}/04-QualityEvaluation
GenomeAnnotation=${ProjectDir}/05-GenomeAnnotation
GenomesComparison=${ProjectDir}/06-GenomesComparison

################################
### Making directory ###
function DirExists(){
	if [ ! -d $1 ]
	then 
		mkdir -p $1
	fi
}

for i in ${Source_data} ${GenomeSurvey} ${Contigs} ${Scaffolds} ${PseudoChromosomes} ${EvaluateQuality} ${GenomeAnnotation}  ${GenomesComparison}
do 
	DirExists $i
done 

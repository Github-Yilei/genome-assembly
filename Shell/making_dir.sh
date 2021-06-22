#usr/bin/bash
echo making dir for genome assembly

################################
### Variable information ###
MainFile=Cclementina
Species=test

################################
### Project information ###
ProjectDir=/share/home/stu_wuyilei/Projects/GenomeAssembly/${Species}

GenomeSurvey=${ProjectDir}/00-GenomeSurvey
Assembly=${ProjectDir}/01-Assembly
EvaluateQuality=${ProjectDir}/02-EvaluateQuality


################################
### Making directory ###
function DirExists(){
	if [ ! -d $1 ]
	then 
		mkdir -p $1
	fi
}

for i in ${RawFastQC} ${CleanFastQC} ${CleanReads} ${FastpDir} ${Log} ${FastpLog} ${PBSLog}
do 
	DirExists $i
done 

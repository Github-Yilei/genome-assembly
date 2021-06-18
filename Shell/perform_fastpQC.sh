#usr/bin/bash
echo start FASTQ quality control and data-filtering with fastp at `date`

################################
### Variable information ###
Species=Species_name
sampleID=reads_prefix

################################
### Project information ###
ProjectDir=/share/home/stu_wuyilei/project/Genome_assembly/${Species}
RawReads=${ProjectDir}/Source_data/Illumina/Raw

RawFastQC=${ProjectDir}/Source_data/Illumina/RawFastQC
CleanFastQC=${ProjectDir}/Source_data/Illumina/CleanFastQC
CleanReads=${ProjectDir}/Source_data/Illumina/Cleaned
FastpDir=${ProjectDir}/Source_data/Illumina/FastpDir

PBSLog=${ProjectDir}/workflow/QsubLog/FastpLog
FastpLog=${Log}/Source_data/Illumina/FastpLog

################################
### Build directory ###
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

echo ${RawFastQC} ...OK
echo ${CleanFastQC} ...OK
echo ${CleanReads} ...OK
echo ${FastpDir} ...OK
echo ${Log} ...OK
echo ${FastpLog} ...OK
echo ${PBSLog} ...OK


################################
### Perform pbs ###
echo everything seems fine, lets submit the Jobs
cd ${ProjectDir}/workflow
echo now, we are in file ${ProjectDir}/workflow, and ready to qsub jobs.

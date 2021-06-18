#usr/bin/bash
echo start FASTQ quality control and data-filtering with fastp at `date`

################################
### Variable information ###
MainFile=MainFile


################################
### Project information ###
ProjectDir=/share/home/stu_wuyilei/project/${MainFile}
SampleList=${ProjectDir}/MetaData/sample_list
RawReads=${ProjectDir}/Rawdata

RawFastQC=${ProjectDir}/FastQCDir/RawFastQC
CleanFastQC=${ProjectDir}/FastQCDir/CleanFastQC
CleanReads=${ProjectDir}/CleanReads

FastpDir=${ProjectDir}/FastpDir

Log=${ProjectDir}/Log
FastpLog=${Log}/FastpLog

PBSLog=${ProjectDir}/workflow/QsubLog/FastpLog
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

### submit missions ###
qsub -J 1-$1:1 fastp_QC.sh

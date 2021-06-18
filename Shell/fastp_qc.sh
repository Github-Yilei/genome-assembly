#!/bin/bash
#PBS -N FastpQC
#PBS -o QsubLog/FastpLog/
#PBS -e QsubLog/FastpLog/
#PBS -q workq
#PBS -l nodes=1:ppn=5
#PBS -A YiLei

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

PBSLog=${ProjectDir}/workflow/QsubLog/FastpLog
FastpLog=${Log}/FastpLog

################################
### Tools ###
Fastp=/share/home/stu_wuyilei/miniconda3/bin/fastp
Fastqc=/share/home/stu_wuyilei/biosoft/FastQC/fastqc
 
################################
### Main parameters ###
RequiredCPU=5

IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###
cd ${ProjectDir}

### FastQC for raw reads ###
${Fastqc} -o ${RawFastQC} -t ${RequiredCPU} ${RawReads}/${sampleID}_1.fq.gz ${RawReads}/${sampleID}_2.fq.gz
echo Fastqc for raw reads finished with ${sampleID} at `date`

### Quality control ###
${Fastp} -i ${RawReads}/${sampleID}_1.fq.gz -I ${RawReads}/${sampleID}_2.fq.gz -o ${CleanReads}/${sampleID}_1.fq.gz \
	-O ${CleanReads}/${sampleID}_2.fq.gz -W 5 -M 20 -5 -3 -l 50 -w ${RequiredCPU} \
	-j ${FastpDir}/${sampleID}.json -h ${FastpDir}/${sampleID}.html > ${FastpLog}/${sampleID}.qc.log 2>&1
echo fastp finished with ${sampleID} at `date`

### Fastqc for clean reads ###
${Fastqc} -o ${CleanFastQC} -t ${RequiredCPU} ${CleanReads}/${sampleID}_1.fq.gz ${CleanReads}/${sampleID}_2.fq.gz
echo Fastqc for clean-reads finished with ${sampleID} at `date`

################################
### adding a finished tag
echo fastp_QC of ${sampleID} is finished at `date`

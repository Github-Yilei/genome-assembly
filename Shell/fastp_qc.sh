#!/bin/bash
#PBS -N FastpQC
#PBS -o QsubLog/FastpLog/
#PBS -e QsubLog/FastpLog/
#PBS -q workq
#PBS -l nodes=1:ppn=30
#PBS -A YiLei

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
### Tools ###
Fastp=/share/home/stu_wuyilei/miniconda3/bin/fastp
Fastqc=/share/home/stu_wuyilei/biosoft/FastQC/fastqc
 
################################
### Main parameters ###
RequiredCPU=30

################################
### Main process ###
cd ${ProjectDir}

### FastQC for raw reads ###
${Fastqc} -o ${RawFastQC} -t ${RequiredCPU} ${RawReads}/${sampleID}_1.fq.gz ${RawReads}/${sampleID}_2.fq.gz
echo Fastqc for raw reads finished with ${sampleID} at `date`

### Quality control ###
${Fastp} -i ${RawReads}/${sampleID}_1.fq.gz -I ${RawReads}/${sampleID}_2.fq.gz -o ${CleanReads}/${sampleID}_1.fq.gz \
	-O ${CleanReads}/${sampleID}_2.fq.gz -W 5 -M 20 -5 -3 -l 50 -z 6 -w ${RequiredCPU} \
	-j ${FastpDir}/${sampleID}.json -h ${FastpDir}/${sampleID}.html > ${FastpLog}/${sampleID}.qc.log 2>&1
echo fastp finished with ${sampleID} at `date`

### Fastqc for clean reads ###
${Fastqc} -o ${CleanFastQC} -t ${RequiredCPU} ${CleanReads}/${sampleID}_1.fq.gz ${CleanReads}/${sampleID}_2.fq.gz
echo Fastqc for clean-reads finished with ${sampleID} at `date`

################################
### adding a finished tag
echo fastp_QC of ${sampleID} is finished at `date`

#!/bin/bash
#PBS -N FastpQC
#PBS -o Fastp.log
#PBS -e Fastp.err
#PBS -q workq
#PBS -l nodes=1:ppn=3
#PBS -A YiLei

################################
### Variable information ###
ProjectDir=Source_data/RNAseq

################################
### Project information ###
RawReads=${ProjectDir}/Data
SampleList=${ProjectDir}/sampleList.txt

RawFastQC=${ProjectDir}/RawFastQC
CleanFastQC=${ProjectDir}/CleanFastQC
CleanReads=${ProjectDir}/Cleaned
FastpDir=${ProjectDir}/FastpDir

FastpLog=${ProjectDir}/FastpLog

################################
### Build directory ###
function DirExists(){
        if [ ! -d $1 ]
        then
                mkdir -p $1
        fi
}

for i in ${RawFastQC} ${CleanFastQC} ${CleanReads} ${FastpDir} ${FastpLog}
do
        DirExists $i
done

echo ${RawFastQC} ...OK
echo ${CleanFastQC} ...OK
echo ${CleanReads} ...OK
echo ${FastpDir} ...OK
echo ${FastpLog} ...OK

################################
### Tools ###
Fastp=/share/home/stu_wuyilei/miniconda3/bin/fastp
Fastqc=/share/home/stu_wuyilei/biosoft/FastQC/fastqc

################################
### Main parameters ###
RequiredCPU=3

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
        -O ${CleanReads}/${sampleID}_2.fq.gz -w ${RequiredCPU} \
        -j ${FastpDir}/${sampleID}.json -h ${FastpDir}/${sampleID}.html > ${FastpLog}/${sampleID}.qc.log 2>&1
echo fastp finished with ${sampleID} at `date`

### Fastqc for clean reads ###
${Fastqc} -o ${CleanFastQC} -t ${RequiredCPU} ${CleanReads}/${sampleID}_1.fq.gz ${CleanReads}/${sampleID}_2.fq.gz
echo Fastqc for clean-reads finished with ${sampleID} at `date`

################################
### adding a finished tag
echo fastp_QC of ${sampleID} is finished at `date`

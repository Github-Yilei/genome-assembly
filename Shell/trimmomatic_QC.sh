#!/bin/bash
#PBS -N TrimmomaticQC
#PBS -o QsubLog/TrimmomaticQC.log
#PBS -e QsubLog/TrimmomaticQC.err
#PBS -q workq
#PBS -l mem=10GB
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

CleanDir=$WorkDir/TrimmomaticClean
TrimmomaticLog=${ProjectDir}/Log/TrimmomaticLog

################################
### Tools ###
Trimmomatic=/share/home/stu_wuyilei/miniconda3/bin/trimmomatic

################################
### Main parameters ###
RequiredCPU=5

IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Build directory ###

function DirExists(){
    if [ ! -d $1 ]
    then
        mkdir -p $1
    fi
}

for i in ${CleanDir} ${TrimmomaticLog}
do
    DirExists $i
done

################################
### Main process ###
cd ${ProjectDir}

### TrimmomaticQC for raw reads ###
echo Trimmomatic Quality Control started with ${sampleID} at `date`

java -jar ${Trimmomatic} PE -threads ${RequiredCPU}  -phred33\
    ${RawReads}/${sampleID}_1.fq.gz \
    ${RawReads}/${sampleID}_2.fq.gz \
    ${CleanDir}/${sampleID}_pair_1.fq.gz \
    ${CleanDir}/${sampleID}_unpair_1.fq.gz \
    ${CleanDir}/${sampleID}_pair_2.fq.gz \
    ${CleanDir}/${sampleID}_unpair_2.fq.gz \
    HEADCROP:15 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:50 \
    -summary ${TrimmomaticLog}/${sampleID}.log

echo Trimmomatic QC finished with ${sampleID} at `date`


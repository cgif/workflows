#!/bin/bash

## script to run GATK

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=8gb

#PBS -m bea
#PBS -j oe

# load modules

module load cnvnator/0.2.5

NOW="date +%Y-%m-%d%t%T%t"

# define variables

SAMPLE=sample
ANALYSIS_DIR=analysisDir
FASTA_FOLDER=fastaFolder
ROOT_FILES="rootFiles"
BIN_SIZE=binSize

# merge root files for all autosomes
echo "`${NOW}`merging root files"

MERGED_ROOT_FILE=${ANALYSIS_DIR}/${SAMPLE}.root

cnvnator -root $MERGED_ROOT_FILE -merge $ROOT_FILES

# generate histogram

cnvnator -root $MERGED_ROOT_FILE -his $BIN_SIZE -d $FASTA_FOLDER

# calculate statistics

cnvnator -root $MERGED_ROOT_FILE -stat $BIN_SIZE

# partition RD signal
# will try for the whole file, if too slow then will do by chromosome

cnvnator -root $MERGED_ROOT_FILE -partition $BIN_SIZE

# cnv calling

cnvnator -root $MERGED_ROOT_FILE -call $BIN_SIZE > ${ANALYSIS_DIR}/${SAMPLE}_cnvs.txt



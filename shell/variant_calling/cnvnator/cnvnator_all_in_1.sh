#!/bin/bash

## script to run GATK

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -m bea
#PBS -j oe

# load modules

module load cnvnator/0.2.5

NOW="date +%Y-%m-%d%t%T%t"

# define variables

INPUT_BAM=inputBam
SAMPLE=sample
ANALYSIS_DIR=analysisDir
FASTA_FOLDER=fastaFolder
BIN_SIZE=binSize

# merge root files for all autosomes
echo "`${NOW}`extracting reads"
ROOT_FILE=$ANALYSIS_DIR/${SAMPLE}.root

cnvnator -root $ROOT_FILE -chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 -tree $INPUT_BAM

# generate histogram
echo "`${NOW}`generating histogram"
cnvnator -root $ROOT_FILE -his $BIN_SIZE -d $FASTA_FOLDER

# calculate statistics
echo "`${NOW}`calculating statistics"
cnvnator -root $ROOT_FILE -stat $BIN_SIZE

# partition RD signal
# will try for the whole file, if too slow then will do by chromosome
echo "`${NOW}`partitioning RD signal"
cnvnator -root $ROOT_FILE -partition $BIN_SIZE

# cnv calling
echo "`${NOW}`calling CNVs"
cnvnator -root $ROOT_FILE -call $BIN_SIZE > ${ANALYSIS_DIR}/${SAMPLE}_cnvs.txt



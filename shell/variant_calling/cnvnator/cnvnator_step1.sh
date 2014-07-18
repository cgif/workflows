#!/bin/bash

## script to run GATK

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=8gb

#PBS -m bea
#PBS -j oe

# load modules

module load cnvnator/0.2.5

NOW="date +%Y-%m-%d%t%T%t"

# define variables

INPUT_BAM=inputBam
SAMPLE=sample
ANALYSIS_DIR=analysisDir
CHR=chromosomeName

#mkdir $TMPDIR/tmp


# extract reads from bam file
echo "`${NOW}`extracting reads from bam file $INPUT_BAM"

ROOT_FILE=${ANALYSIS_DIR}/per_chromosome/${SAMPLE}.chr${CHR}.root

cnvnator -root $ROOT_FILE -chrom $CHR -tree $INPUT_BAM


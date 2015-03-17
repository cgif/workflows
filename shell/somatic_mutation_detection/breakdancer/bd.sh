#!/bin/bash

## script to run breakDancer for translocations

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb:tmpspace=#tmpSpacegb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

# load modules
module load samtools/#samtoolsVersion
module load breakdancer/#bdVersion

#now
NOW="date +%Y-%m-%d%t%T%t"

# define variables
INPUT_BAM=#inputBam
RESULTS_FOLDER=#resultsFolder

BAM_NAME=`basename $INPUT_BAM .bam`

#copy bam files into TMPDIR
cp $INPUT_BAM $TMPDIR/${BAM_NAME}.bam
cp $INPUT_BAM.bai $TMPDIR/${BAM_NAME}.bam.bai

#create configuration file
echo "`${NOW}`create configuration file"
bam2cfg.pl -g $TMPDIR/${BAM_NAME}.bam > $TMPDIR/${BAM_NAME}.cfg

#detect SVs
echo "`${NOW}`detect SVs"
breakdancer_max -t -r 3 -g $TMPDIR/${BAM_NAME}.bed $TMPDIR/${BAM_NAME}.cfg > $TMPDIR/${BAM_NAME}.bd

#copy output files to results folder
cp $TMPDIR/${BAM_NAME}.bed $RESULTS_FOLDER
cp $TMPDIR/${BAM_NAME}.bd $RESULTS_FOLDER

chmod 0770 $RESULTS_FOLDER
chmod 0660 $RESULTS_FOLDER/*

echo "`${NOW}` done"


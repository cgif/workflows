#!/bin/bash

## script to run CREST

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

# load modules
module load samtools/0.1.18
module load breakdancer

#now
NOW="date +%Y-%m-%d%t%T%t"

# define variables

INPUT_BAM=#inputBam
RESULTS_FOLDER=#resultsFolder

BAM_NAME=`basename $INPUT_BAM .bam`

#copy bam files into TMPDIR
cp $INPUT_BAM $TMPDIR/${BAM_NAME}.bam

#create configuration file
echo "`${NOW}`create configuration file"
/apps/breakdancer/2012-07-04/perl/bam2cfg.pl -g -h $TMPDIR/${BAM_NAME}.bam > $TMPDIR/${BAM_NAME}.cfg

#detect SVs
echo "`${NOW}`detect SVs"
/apps/breakdancer/2012-07-04/breakdancer -g $TMPDIR/${BAM_NAME}.bed $TMPDIR/${BAM_NAME}.cfg > $TMPDIR/${BAM_NAME}.bd

#copy output files to results folder
cp $TMPDIR/${BAM_NAME}.cfg $RESULTS_FOLDER
cp $TMPDIR/${BAM_NAME}.bed $RESULTS_FOLDER
cp $TMPDIR/${BAM_NAME}.bd $RESULTS_FOLDER

chmod 770 $RESULTS_FOLDER
chmod 660 $RESULTS_FOLDER/*

echo "`${NOW}` done"


#!/bin/bash

#run R script to create MEDIPS sets and
#find regions of differential methylation

#PBS -l walltime=72:00:00
#PBS -l ncpus=1
#PBS -l mem=100gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/#Rversion

R_SCRIPT="#Rscript"
INPUT_DIR="#inputBam"
SAMPLE_INFO="#sampleInfo"
CONDITION1="#condition1"
CONDITION2="#condition2"

#copy bam files into tmp space
for SAMPLE in `sed 1d $SAMPLE_INFO | grep -Pi "$CONDITION1|$CONDITION2|input" | cut -f1 | grep -vP "^$"`;do

	echo "`${NOW}`copying bam file for sample $SAMPLE"
	BAM_PATH=$INPUT_DIR/$SAMPLE/$SAMPLE.nondup.rename.filt.bam
	BAM_TMP=$TMPDIR/$SAMPLE.bam
	cp $BAM_PATH $BAM_TMP

done

sed -i -e "s/#tmpDir/${TMPDIR//\//\\/}/" $R_SCRIPT

#calculate methylation profiles
R --vanilla < $R_SCRIPT

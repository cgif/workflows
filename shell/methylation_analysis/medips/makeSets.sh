#!/bin/bash

#run R script to create MEDIPS sets,
#calculate methylation profiles, 
#make WIG files to open in UCSC browser

#PBS -l walltime=72:00:00
#PBS -l ncpus=1
#PBS -l mem=100gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/#Rversion

R_SCRIPT="#Rscript"
INPUT_DIR="#inputBam"
SAMPLES="#samples"
CONDITIONS="#conditions"

#copy bam files into tmp space
for SAMPLE in $SAMPLES; do

	echo "`${NOW}`copying bam file for sample $SAMPLE"
	BAM_PATH=$INPUT_DIR/$SAMPLE/$SAMPLE.nondup.rename.filt.bam
	BAM_TMP=$TMPDIR/$SAMPLE.bam
	cp $BAM_PATH $BAM_TMP

done

sed -i -e "s/#tmpDir/${TMPDIR//\//\\/}/" $R_SCRIPT

#calculate methylation profiles
R --vanilla < $R_SCRIPT

#make WIG files from methylation profiles
#for CONDITION in $CONDITIONS; do

#	WIG=""
#	PROFILE=""
#	echo "track type=wiggle_0 name=$condition.$set_up description=$condition.$set_up visibility=full autoScale=on color=0,0,255 maxHeightPixels=100:50:20 graphType=bar priority=20\n" > $WIG

#	for CHROM in `cut -f 2 $PROFILE|sort|uniq`; do

#		echo "fixedStep chrom=$CHROM start=1 step=100 span=100\n" > $WIG
#		cat $WIG|awk ($2 == $CHROM) print|cut -f 10 > $WIG

#	done
#done

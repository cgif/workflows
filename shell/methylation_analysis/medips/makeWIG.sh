#!/bin/bash

#run R script to create MEDIPS sets,
#calculate methylation profiles,
#make bigWig files and open them in UCSC browser

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
MEDIPS_DIR="#medipsDir"
CHROM_SIZES="#chromSizes"
CONDITION="#condition"
WINDOW=#windowSize
ENC_DIR="#encDir"
DEPLOYMENT_SERVER="#deploymentServer"

#copy bam files into tmp space
for SAMPLE in `sed 1d $SAMPLE_INFO | grep "$CONDITION" | cut -f1 | grep -vP "^$"`;do

	echo "`${NOW}`copying bam file for sample $SAMPLE"
	BAM_PATH=$INPUT_DIR/$SAMPLE/$SAMPLE.hg19.bam
	BAM_TMP=$TMPDIR/$SAMPLE.bam
	cp $BAM_PATH $BAM_TMP

done

sed -i -e "s/#tmpDir/${TMPDIR//\//\\/}/" $R_SCRIPT

#calculate methylation profiles
R --vanilla < $R_SCRIPT

#make bigWig files from methylation profiles
PROFILE="$MEDIPS_DIR/WIG/${CONDITION}_profile.txt"
BIG_WIG="$MEDIPS_DIR/WIG/${CONDITION}.bw"

TMP_PROFILE="$TMPDIR/${CONDITION}_profile.txt"
TMP_WIG="$TMPDIR/${CONDITION}.wig"
TMP_BIG_WIG="$TMPDIR/${CONDITION}.bw"

cp $PROFILE $TMP_PROFILE
echo -n "" > $TMP_WIG

#create WIG file from the last column of methylation profile
for CHROM in `cut -f 1 $CHROM_SIZES`; do

	echo "fixedStep chrom=$CHROM start=1 step=$WINDOW span=$WINDOW" >> $TMP_WIG
	cat $TMP_PROFILE|sed -e 's/"//g'|awk -v CHROM=$CHROM '{if ($2 == CHROM) print}'|rev|cut -f 1|rev >> $TMP_WIG

done

/groupvol/cgi/software/ucsc/wigToBigWig $TMP_WIG $CHROM_SIZES $TMP_BIG_WIG
cp $TMP_BIG_WIG $BIG_WIG
scp $BIG_WIG $DEPLOYMENT_SERVER:$ENC_DIR

	
	







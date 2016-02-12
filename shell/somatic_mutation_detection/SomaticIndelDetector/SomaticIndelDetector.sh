#!/bin/bash

## script to run SomaticIndelDetector

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules
module load java/#javaVersion
module load gatk/2.3-9

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

# define variables
ANALYSIS_FILE=#analysisFile
REFERENCE_FASTA=#referenceFasta
REFERENCE_DICT=#referenceDict
NORMAL_BAM=#normalBam
TUMOR_BAM=#tumorBam
INTERVALS_FILE=#intervalsFile
SUMMARY_SCRIPT=#summaryScriptPath

#copy input files to tmp dir
echo "`${NOW}` copying files to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/ref.fasta
cp $REFERENCE_FASTA.fai $TMPDIR/ref.fasta.fai
cp $REFERENCE_DICT $TMPDIR/ref.dict
cp $NORMAL_BAM $TMPDIR/normal.bam
cp $NORMAL_BAM.bai $TMPDIR/normal.bam.bai
cp $TUMOR_BAM $TMPDIR/tumor.bam
cp $TUMOR_BAM.bai $TMPDIR/tumor.bam.bai
cp $INTERVALS_FILE $TMPDIR/intervals.intervals

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

echo "`${NOW}` running  SomaticIndelDetector..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
-T SomaticIndelDetector \
-I:normal $TMPDIR/normal.bam \
-I:tumor $TMPDIR/tumor.bam \
-R $TMPDIR/ref.fasta \
-L $TMPDIR/intervals.intervals \
-o $TMPDIR/tmp.vcf \
-filter 'T_COV<14||N_COV<8||T_INDEL_F<0.15||T_INDEL_CF<0.5' \
-ws 200 \
-verbose $TMPDIR/tmp.stats

echo "`${NOW}` copying files from tmp directory..."
cp $TMPDIR/tmp.vcf $ANALYSIS_FILE.vcf
cp $TMPDIR/tmp.stats $ANALYSIS_FILE.stats

perl $SUMMARY_SCRIPT
#!/bin/bash

## script to run mutect

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules
module load mutect/#mutectVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

# define variables
ANALYSIS_FILE=#analysisFile
REFERENCE_FASTA=#referenceFasta
REFERENCE_DICT=#referenceDict
MUTECT_COSMIC=#mutectCosmic
MUTECT_DBSNP=#mutectDBsnp
NORMAL_BAM=#normalBam
TUMOR_BAM=#tumorBam
INTERVALS_FILE=#intervalsFile
SUMMARY_SCRIPT_PATH=#summaryScriptPath

#copy input files to tmp dir
echo "`${NOW}` copying files to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/ref.fasta
cp $REFERENCE_FASTA.fai $TMPDIR/ref.fasta.fai
cp $REFERENCE_DICT $TMPDIR/ref.dict
cp $MUTECT_COSMIC $TMPDIR/cosmic.vcf
cp $MUTECT_DBSNP $TMPDIR/dbsnp.vcf
cp $NORMAL_BAM $TMPDIR/normal.bam
cp $NORMAL_BAM.bai $TMPDIR/normal.bam.bai
cp $TUMOR_BAM $TMPDIR/tumor.bam
cp $TUMOR_BAM.bai $TMPDIR/tumor.bam.bai
cp $INTERVALS_FILE $TMPDIR/intervals.intervals

# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

echo "`${NOW}` running mutect..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -jar -Djava.io.tmpdir=$TMPDIR/tmp $MUTECT_HOME/muTect-1.1.4.jar \
--analysis_type MuTect \
--reference_sequence $TMPDIR/ref.fasta \
--cosmic $TMPDIR/cosmic.vcf \
--dbsnp $TMPDIR/dbsnp.vcf \
--intervals $TMPDIR/intervals.intervals \
--input_file:normal $TMPDIR/normal.bam \
--input_file:tumor $TMPDIR/tumor.bam \
--out $TMPDIR/tmp.txt \
--vcf $TMPDIR/tmp.vcf \
-rf BadCigar

echo "`${NOW}` copying files from tmp directory..."
cp $TMPDIR/tmp.txt $ANALYSIS_FILE.stats
cp $TMPDIR/tmp.vcf $ANALYSIS_FILE.vcf

perl $SUMMARY_SCRIPT_PATH
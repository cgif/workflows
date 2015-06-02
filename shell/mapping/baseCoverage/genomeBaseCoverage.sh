#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load gatk/3.3
module load java/jdk-7u25

JAVA_XMX=4G

# define variables
REFERENCE_FASTA=#pathReferenceFasta
INPUT_BAM=#inputBam
RESULTS_DIR=#resultsDir
SAMPLE=#sampleName

mkdir $TMPDIR/tmp

echo "`${NOW}`INFO copying input BAM and index files to tmp directory..."
cp $INPUT_BAM $TMPDIR/input.bam
cp $INPUT_BAM.bai $TMPDIR/input.bam.bai 

echo "`${NOW}`INFO $SCRIPT_CODE copying reference to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict


#########################################################
# run BaseCoverageDistribution
#########################################################

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-T BaseCoverageDistribution \
	-R $TMPDIR/reference.fa \
	-I $TMPDIR/input.bam \
	-rf BadCigar \
	-o $TMPDIR/report.grp


echo "`${NOW}`INFO copying report.grp to output directory $RESULTS_DIR..."
cp $TMPDIR/report.grp $RESULTS_DIR/$SAMPLE.report.grp

ls -al $TMPDIR

echo "`${NOW}`INFO done"


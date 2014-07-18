#!/bin/bash

## test script to prepare bam for pindel analysis
## removes duplicates from bam and collects insert metrics

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=6gb

#PBS -m bea
#PBS -j oe

# load modules

module load samtools/0.1.18
module load picard/1.65
module load R/2.15

NOW="date +%Y-%m-%d%t%T%t"
JAVA_MXM=5800M

# define variables

INPUT_BAM=inputBam
ANALYSIS_DIR=analysisDir
RESULTS_DIR=resultsDir
SAMPLE=sample

mkdir $TMPDIR/tmp

# step 1: remove duplicates from bam file
echo "`${NOW}`removing flagged duplicates"
samtools view -b -F 0x400 -o $ANALYSIS_DIR/${SAMPLE}.recalibrated.unique.bam $INPUT_BAM 

# step 2: index dedupped bam file
echo "`${NOW}`indexing dedupped bam file"
samtools index $ANALYSIS_DIR/${SAMPLE}.recalibrated.unique.bam

# step 3: calculate isert size metrics using picard
echo "`${NOW}`calculating insert size metrics"
java -Xmx$JAVA_MXM -jar $PICARD_HOME/CollectInsertSizeMetrics.jar \
  TMP_DIR=$TMPDIR/tmp \
  MAX_RECORDS_IN_RAM=1000000 \
  INPUT=$ANALYSIS_DIR/${SAMPLE}.recalibrated.unique.bam \
  HISTOGRAM_FILE=$ANALYSIS_DIR/${SAMPLE}.isize.histogram.pdf \
  METRIC_ACCUMULATION_LEVEL=SAMPLE \
  OUTPUT=$ANALYSIS_DIR/${SAMPLE}.isizeMetrics \
  VALIDATION_STRINGENCY=SILENT

cp $ANALYSIS_DIR/${SAMPLE}.isizeMetrics $RESULTS_DIR
cp $ANALYSIS_DIR/${SAMPLE}.isize.histogram.pdf $RESULTS_DIR

echo "`${NOW}`done"

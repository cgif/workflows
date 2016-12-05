#!/bin/bash

## script to bismark alignment and methylation extractor
## written for a small dataset which can be processed all at once
## for large datasets the scripts may need to be rewritten 
## to use more cores and to split/merge files
## bismark by default uses several cores (at least 5 with bowtie2)

#PBS -l walltime=#walltimeHours:00:00
#PBS -l select=1:ncpus=#threads:mem=50gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

# load modules

module load bismark/#bismarkVersion
module load bowtie/#bowtieVersion
module load samtools/#samtoolsVersion

NOW="date +%Y-%m-%d%t%T%t"

REFERENCE_DIR=#referenceDir
RESULTS_DIR=#pathOutputDir

PATH_READS_DIRECTORY=#pathReadsDirectory
FASTQ_READ1=#read1
FASTQ_READ2=#read2
SAMPLE=#sampleName
CUSTOMER_ID=#customerID

echo "`${NOW}`copying references to tmp directory..."
cp -R $REFERENCE_DIR $TMPDIR

echo "`${NOW}`copying data files to tmp directory..."

cp $PATH_READS_DIRECTORY/$FASTQ_READ1 $TMPDIR/R1.fq.gz
cp $PATH_READS_DIRECTORY/$FASTQ_READ2 $TMPDIR/R2.fq.gz

#make temporary output directories
mkdir $TMPDIR/bam
mkdir $TMPDIR/methylation

echo "`${NOW}`running bismark alignment"

bismark --genome $TMPDIR/fasta \
	-1 $TMPDIR/R1.fq.gz \
	-2 $TMPDIR/R2.fq.gz \
	--basename $SAMPLE \
	--rg_tag TRUE \
	--rg_id $SAMPLE \
	--rg_sample $CUSTOMER_ID \
	-o $TMPDIR/bam

echo "`${NOW}`sort and index bam file"
samtools sort $TMPDIR/bam/${SAMPLE}_pe.bam $TMPDIR/bam/${SAMPLE}_pe.sorted
mv $TMPDIR/bam/${SAMPLE}_pe.sorted.bam $TMPDIR/bam/${SAMPLE}_pe.bam
samtools index $TMPDIR/bam/${SAMPLE}_pe.bam

echo "`${NOW}`contents of alignment directory on $TMPDIR/bam"
ls -al $TMPDIR/bam

echo "`${NOW}`copying alignment to results directory"
cp -R $TMPDIR/bam $RESULTS_DIR

#run methylation extractor

bismark_methylation_extractor $TMPDIR/bam/${SAMPLE}_pe.bam \
	-o $TMPDIR/methylation \
	--bedGraph

echo "`${NOW}`contents of methylation directory on $TMPDIR/methylation"
ls -al $TMPDIR/methylation

echo "`${NOW}`copying output to results directory"
cp -R $TMPDIR/methylation $RESULTS_DIR
echo "`${NOW}`done"

echo "`${NOW}`contents of input directory on $TMPDIR"
ls -al $TMPDIR

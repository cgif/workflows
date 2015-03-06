#!/bin/bash

#
# script to generate cram files from bam files
# configured by submitMapping script
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=threads:mem=7800mb:tmpspace=50gb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi




module load samtools/1.1

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#bwa aln number of threads
THREADS=threads

#input BAM
PATH_INPUT_BAM=pathInputBam

#cram output directory
PATH_OUTPUT_CRAM_DIR=pathOutputCramDir

#output prefix
OUTPUT_PREFIX=outputPrefix

						
#CRAM_NAME=`echo $PATH_INPUT_BAM | awk 'BEGIN{FS="/"} {print $8}' | awk 'BEGIN {FS="."} {print $1}'`
CRAM_NAME=`echo $OUTPUT_PREFIX | awk 'BEGIN {FS="."} {print $1}'`
PATH_OUTPUT_CRAM=$PATH_OUTPUT_CRAM_DIR/$CRAM_NAME.cram

#path to reference genome fasta file without gzip extension
PATH_REFERENCE_FASTA_NO_EXT=pathReferenceFastaNoExt
REFERENCE_FASTA_NAME=`basename $PATH_REFERENCE_FASTA_NO_EXT`

echo "`${NOW}`staging input files..." 
echo "`${NOW}`$PATH_INPUT_BAM" 
cp $PATH_INPUT_BAM $TMPDIR/tmp.bam  
echo "`${NOW}`$PATH_REFERENCE_FASTA_NO_EXT"
cp $PATH_REFERENCE_FASTA_NO_EXT $TMPDIR/$REFERENCE_FASTA_NAME

echo "`${NOW}`converting BAM to CRAM..." 
samtools view -@ $THREADS -T $TMPDIR/$REFERENCE_FASTA_NAME -C -o $TMPDIR/tmp.cram $TMPDIR/tmp.bam

echo "`${NOW}`done" 

echo "`${NOW}`listing the files in temporary folder for debugging"
echo "`ls -l $TMPDIR`"

echo "`${NOW}`copying CRAM to $PATH_OUTPUT_CRAM" 
cp $TMPDIR/tmp.cram $PATH_OUTPUT_CRAM
chmod 660 $PATH_OUTPUT_CRAM


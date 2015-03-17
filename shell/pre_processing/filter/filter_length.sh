#!/bin/bash

#
# script filter out short sequences after cutadapt
#

#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:mem=4gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

#path to fastq file 

PATH_READS_FASTQ=pathReadsFastq

MIN_LENGTH=#minLength


FASTQ_NAME=`basename $PATH_READS_FASTQ`
FASTQ_BASENAME=`echo $FASTQ_NAME | perl -pe 's/(.*)\.f.*q.*/$1/g'`

#copy files to temp directory
echo "`${NOW}`copying read $FASTQ_NAME to temporary scratch space..."
cp $PATH_READS_FASTQ $TMPDIR/$FASTQ_NAME

echo "`${NOW}`listing $TMPDIR content"
ls $TMPDIR

echo "`${NOW}`unzipping $TMPDIR/$FASTQ_NAME"
FASTQ_UNZIPPED=`basename $TMPDIR/$FASTQ_NAME .gz`
gunzip $TMPDIR/$FASTQ_NAME
echo "`${NOW}`listing $TMPDIR content"
ls $TMPDIR

awk 'BEGIN {OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= #minLength ) {print header, seq, qheader, qseq}}' < $TMPDIR/$FASTQ_UNZIPPED > $TMPDIR/${FASTQ_BASENAME}_filtered.fq

gzip $TMPDIR/${FASTQ_BASENAME}_filtered.fq

OUTPUT_PATH=`dirname $PATH_READS_FASTQ`

echo "`${NOW}`copying filtered fastqc to $OUTPUT_PATH"
cp $TMPDIR/${FASTQ_BASENAME}_filtered.fq.gz $OUTPUT_PATH
chmod 660 $OUTPUT_PATH/*

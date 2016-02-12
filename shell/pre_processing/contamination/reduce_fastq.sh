#!/bin/bash

#
# script to extract fraction of reads from fastq file
#

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

module load python/#python_version

COUNT=#count
FASTQ1=#fastq1_path
FASTQ2=#fastq2_path
OUTPUT_DIR=#analysis_dir
BASEDIR=#base_dir

LINE_COUNT=`gzip -dc $FASTQ1 | wc -l`
TOTAL_READS=$((LINE_COUNT/4))
FRACTION=`echo "scale=2;$COUNT/$TOTAL_READS"|bc`

FASTQ1_BASENAME=`basename $FASTQ1`
FASTQ2_BASENAME=`basename $FASTQ2`

gzip -c -d $FASTQ1 > $TMPDIR/fastq1.fq
gzip -c -d $FASTQ2 > $TMPDIR/fastq2.fq

python $BASEDIR/reduce_fastq.py $FRACTION $TMPDIR/fastq1.fq $TMPDIR/fastq2.fq $TMPDIR/fastq1_reduced.fq $TMPDIR/fastq2_reduced.fq

gzip -c $TMPDIR/fastq1_reduced.fq > $OUTPUT_DIR/$FASTQ1_BASENAME
gzip -c $TMPDIR/fastq2_reduced.fq > $OUTPUT_DIR/$FASTQ2_BASENAME

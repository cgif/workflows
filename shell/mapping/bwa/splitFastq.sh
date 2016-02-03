#!/bin/bash

#
# script to split fastq files
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l select=1:ncpus=#threads:mem=1024mb:tmpspace=#tmpSpacemb
#PBS -l place=excl

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi


#now
NOW="date +%Y-%m-%d%t%T%t"

#parallel gzip home
PIGZ_HOME=/project/tgu/software/pigz-2.3.1

INPUT_FASTQ=#inputFastq
INPUT_FASTQ_NAME=`basename $INPUT_FASTQ .gz`
READS_PER_CHUNK=#readsPerChunk
LINES_PER_FILE=$(($READS_PER_CHUNK * 4))
OUTPUT_DIR=#outputDir
PIGZ_THREADS=#threads

#copy input fastq file to node
echo "`$NOW`copying fastq file to $TMPDIR... "
cp $INPUT_FASTQ $TMPDIR/

#create temporary directory for split files
mkdir -p $TMPDIR/split

#uncompress and split fastq file
#$PIGZ_HOME/pigz -p $PIGZ_THREADS -cd $INPUT_FASTQ_NAME \
#		| split -d -l $LINES_PER_FILE - \
#		$TMPDIR/split/$INPUT_FASTQ_NAME.

echo "`$NOW`uncompressing fastq file... "	   	
$PIGZ_HOME/pigz -p $PIGZ_THREADS -d $INPUT_FASTQ_NAME.gz \

echo "`$NOW`splitting fastq files into chunks of $READS_PER_CHUNK reads... "
split -d -l $LINES_PER_FILE $INPUT_FASTQ_NAME \
		$TMPDIR/split/$INPUT_FASTQ_NAME.

echo "`$NOW`copying split fastq files to $OUTPUT_DIR... "

if [ -d $OUTPUT_DIR ]; then

	cp $TMPDIR/split/* $OUTPUT_DIR/

else

	mkdir -p $OUTPUT_DIR
	cp $TMPDIR/split/* $OUTPUT_DIR/

fi

ls -al $TMPDIR/
ls -al $TMPDIR/split


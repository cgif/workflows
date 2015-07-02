#!/bin/bash

# script run cutadapt algorithm for removing primer sequences from the ends of the reads

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=4gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load python/2.7.3    # this has cutadapt 1.7.1

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#define variables
#path to fastq file 
PATH_READS_FASTQ=pathReadsFastq

#path to trimmed fastq file
PATH_TRIMMED_DIR=pathTrimmedDir

#primer sequences
PRIMERS_FILE=primersFile

#mkdir -p $PATH_TRIMMED_DIR

#cutadapt parameters

FASTQ_NAME=`basename $PATH_READS_FASTQ`
FASTQ_BASENAME=`echo $FASTQ_NAME | perl -pe 's/(.*)\.f.*q.*/$1/g'`

#copy files to temp directory
echo "`${NOW}`copying read $FASTQ_NAME to temporary scratch space..."
cp $PATH_READS_FASTQ $TMPDIR/$FASTQ_NAME
cp $PRIMERS_FILE $TMPDIR/primers.txt

# remove empty lines from primer file
perl -pi -e 's/^\s*\n//' primers.txt

N_PRIMERS=`cat primers.txt | sed '/^\s*$/d' | wc -l`   

COUNT=1
#rename input to add counter
mv $TMPDIR/$FASTQ_NAME $TMPDIR/$FASTQ_BASENAME.$COUNT.fq.gz

while read LINE; do
#	echo "count at the start $COUNT"
	PRIMER="-g ^$LINE"

#	echo "$PRIMER"
	INPUT_COUNT=$COUNT
	OUTPUT_COUNT=`expr $COUNT + 1`

	echo "`${NOW}`running cutadapt for primer no $OUTPUT_COUNT $LINE"
	cutadapt $PRIMER -o $TMPDIR/$FASTQ_BASENAME.$OUTPUT_COUNT.fq.gz $TMPDIR/$FASTQ_BASENAME.$INPUT_COUNT.fq.gz

	if [[ "$COUNT" -lt "$N_PRIMERS" ]]; then
		echo "deleting intermediate input file $TMPDIR/$FASTQ_BASENAME.$INPUT_COUNT.fq.gz"
		rm $TMPDIR/$FASTQ_BASENAME.$INPUT_COUNT.fq.gz
	else
		echo "$TMPDIR/$FASTQ_BASENAME.$OUTPUT_COUNT.fq.gz is final output"
	fi

	COUNT=$(($COUNT+1))
#	echo "count at the end $COUNT"
done < $TMPDIR/primers.txt

echo "`${NOW}`copying trimmed fastqc $FASTQ_BASENAME.$COUNT.fq.gz to $PATH_TRIMMED_DIR"
cp $TMPDIR/$FASTQ_BASENAME.$COUNT.fq.gz $PATH_TRIMMED_DIR/${FASTQ_BASENAME}_trimmed.fq.gz
chmod 660 $PATH_TRIMMED_DIR/*



#!/bin/bash

#
# script to run BWA alignment for paired end short reads
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l select=1:ncpus=#threads:mem=7800mb:tmpspace=50gb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load bio-bwa/0.7.10
module load samtools/1.1

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"

#bwa aln number of threads
THREADS=#threads

#bwa data directory
PATH_OUTPUT_DIR=#pathOutputDir

#path to reference genome fasta file without gzip extension
PATH_REFERENCE_FASTA_NO_EXT=#pathReferenceFastaNoExt

#path reference fasta name
REFERENCE_FASTA_NAME=`basename $PATH_REFERENCE_FASTA_NO_EXT`

#path to reference genome sequence dictionary
PATH_REFERENCE_DICT=#pathReferenceDict

#path to reference index
PATH_REFERENCE_INDEX_DIR=#pathReferenceIdxDir

#path to reads fastq file without gzip extension
PATH_READS_FASTQ1_NO_EXT=#pathReadsFastqRead1NoExt
PATH_READS_FASTQ2_NO_EXT=#pathReadsFastqRead2NoExt

#output prefix
OUTPUT_PREFIX=#outputPrefix

#keep multiple alignment hits
MULT_READS=#multReads

#path to alignment output file
 
PATH_ALN_BAM=$PATH_OUTPUT_DIR/$OUTPUT_PREFIX.unsorted.bam

TMP_PATH_REFERENCE_FASTA_NO_EXT=$TMPDIR/tmp_reference.fa

TMP_PATH_ALN_BAM_PREFIX=$TMPDIR/tmp_aln.unsorted

#copy reference sequence data
#fasta file is not needed for bwa as it uses the pac file
echo "`${NOW}`copying reference to temporary scratch space..."
for EXT in amb ann bwt pac sa
do
  cp $PATH_REFERENCE_INDEX_DIR/$REFERENCE_FASTA_NAME.$EXT $TMP_PATH_REFERENCE_FASTA_NO_EXT.$EXT
done;

#copy read data to temporary file paths
#paired templorary fastq file is kept empty if data is single-ended
echo "`${NOW}`copying reads to temporary scratch space..."
TMP_PATH_READS_FASTQ1_SUBSET=$TMPDIR/tmp_reads1_subset.fq
cp $PATH_READS_FASTQ1_NO_EXT $TMP_PATH_READS_FASTQ1_SUBSET

if [ -s $PATH_READS_FASTQ2_NO_EXT ]; then
	TMP_PATH_READS_FASTQ2_SUBSET=$TMPDIR/tmp_reads2_subset.fq
	cp $PATH_READS_FASTQ2_NO_EXT $TMP_PATH_READS_FASTQ2_SUBSET
else
	TMP_PATH_READS_FASTQ2_SUBSET=""
fi

#count input reads
echo -n "`${NOW}`counting reads (read1 + read2) in input files: "
INPUT_LINE_COUNT=`wc -l $TMP_PATH_READS_FASTQ1_SUBSET | cut -f1 -d ' '`
#devide line number by 2 (not by 4 because we want the total
#number of reads (read1+read2) not the total number of read
#pairs
#devide by 2 if data is single-ended
if [ -s $PATH_READS_FASTQ2_NO_EXT ]; then
	INPUT_READ_COUNT=$(($INPUT_LINE_COUNT/2))
else
	INPUT_READ_COUNT=$(($INPUT_LINE_COUNT/4))
fi
echo $INPUT_READ_COUNT


###################
#run BWA alignment 
###################

echo "`${NOW}`running BWA alignment"

if [ "$MULT_READS" = "T" ]; then
    echo "`${NOW}` bwa produce multiple hits"
    BWA_CMD="bwa mem -M -t $THREADS $TMP_PATH_REFERENCE_FASTA_NO_EXT $TMP_PATH_READS_FASTQ1_SUBSET $TMP_PATH_READS_FASTQ2_SUBSET | samtools view -bS - > $TMP_PATH_ALN_BAM_PREFIX.bam"
    bwa mem -M -t $THREADS $TMP_PATH_REFERENCE_FASTA_NO_EXT $TMP_PATH_READS_FASTQ1_SUBSET $TMP_PATH_READS_FASTQ2_SUBSET | samtools view -bS - > $TMP_PATH_ALN_BAM_PREFIX.bam				#CHANGE 6 -added this

else

    BWA_CMD="bwa mem -M -t $THREADS $TMP_PATH_REFERENCE_FASTA_NO_EXT $TMP_PATH_READS_FASTQ1_SUBSET $TMP_PATH_READS_FASTQ2_SUBSET | samtools view -bS -F 256 - > $TMP_PATH_ALN_BAM_PREFIX.bam"
    bwa mem -M -t $THREADS $TMP_PATH_REFERENCE_FASTA_NO_EXT $TMP_PATH_READS_FASTQ1_SUBSET $TMP_PATH_READS_FASTQ2_SUBSET | samtools view -bS -F 256 - > $TMP_PATH_ALN_BAM_PREFIX.bam			#CHANGE 6 -added this

fi

#execute command
#$BWA_CMD

#make SAM header

#SAM format version
if [ -e $PATH_REFERENCE_DICT ]
then
	
	#get @HD line of SAM header form BAM file
	samtools view -H $TMP_PATH_ALN_BAM_PREFIX.bam | grep '@HD' > $TMPDIR/header.sam

	#append SN entries from reference dictionary (skipping first line of dictionary)
	cat $PATH_REFERENCE_DICT | grep '@SQ' >> $TMPDIR/header.sam

else

	samtools view -H $TMP_PATH_ALN_BAM_PREFIX.bam > $TMPDIR/header.sam

fi

#program information
echo -e "@PG\tID:BWA\tPN:mem\tVN:$BWA_VERSION_ALIGN\tCL:$BWA_CMD" >> $TMPDIR/header.sam

#replace header
samtools reheader $TMPDIR/header.sam $TMP_PATH_ALN_BAM_PREFIX.bam > $TMP_PATH_ALN_BAM_PREFIX.bam.tmp 
mv $TMP_PATH_ALN_BAM_PREFIX.bam.tmp $TMP_PATH_ALN_BAM_PREFIX.bam



echo "`${NOW}`listing the files in temporary folder for debugging"
echo "`ls -l $TMPDIR`"

#make sure that output has correct number of reads
OUTPUT_READ_COUNT=`samtools flagstat $TMP_PATH_ALN_BAM_PREFIX.bam | head -n 1 | cut -f 1 -d ' '`

echo "`${NOW}`input reads (read1 + read2): $INPUT_READ_COUNT"
echo "`${NOW}`reads in output: $OUTPUT_READ_COUNT"

if [ $OUTPUT_READ_COUNT -ge $INPUT_READ_COUNT ]
then

	echo "`${NOW}`deleting temporary fastq files of read subset..."
	echo "`${NOW}`$PATH_READS_FASTQ1_NO_EXT"	
	rm $PATH_READS_FASTQ1_NO_EXT
	if [ -s $PATH_READS_FASTQ2_NO_EXT ]; then
		echo "`${NOW}`$PATH_READS_FASTQ2_NO_EXT"
        	rm $PATH_READS_FASTQ2_NO_EXT
	fi

else

	echo "`${NOW}`WARNING!!! Output BAM contains less reads then the input fastq file!"
	echo "`${NOW}`keeping temporary fastq files of read subset for re-run" 	

fi

echo "`${NOW}`copying ${TMP_PATH_ALN_BAM_PREFIX}.bam to $PATH_ALN_BAM"
cp $TMP_PATH_ALN_BAM_PREFIX.bam $PATH_ALN_BAM
chmod 640 $PATH_ALN_BAM


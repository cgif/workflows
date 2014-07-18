#!/bin/bash

#
# script to run TopHat alignment for paired end short reads
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l select=1:ncpus=#threads:mem=50gb:tmpspace=25gb
#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

#load modules
module load tophat/#tophatVersion
module load bowtie/#bowtieVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion
module load picard/#picardVersion

JAVA_XMX=19g

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#number of threads
THREADS=#threads

#path to reference genome fasta file without gzip extension
PATH_REFERENCE_FASTA_NO_EXT=#pathReferenceFastaNoExt
REFERENCE_FASTA_BASENAME=`basename $PATH_REFERENCE_FASTA_NO_EXT`

echo "`${NOW}`copying reference and indices to temporary scratch space..."
cp $PATH_REFERENCE_FASTA_NO_EXT.fa $TMPDIR/
cp $PATH_REFERENCE_FASTA_NO_EXT.*.bt2 $TMPDIR/
cp $PATH_REFERENCE_FASTA_NO_EXT.dict $TMPDIR/

#path to annotation file index
PATH_ANNOTATION_GFF=#pathAnnotation
PATH_ANNOTATION_NO_EXT=${PATH_ANNOTATION_GFF%.*}
ANNOTATION_BASENAME=`basename $PATH_ANNOTATION_NO_EXT`

echo "`${NOW}`copying annotation to temporary scratch space..."
cp $PATH_ANNOTATION_NO_EXT.gff $TMPDIR/
cp $PATH_ANNOTATION_NO_EXT.fa* $TMPDIR/
cp $PATH_ANNOTATION_NO_EXT.ver $TMPDIR/
cp $PATH_ANNOTATION_NO_EXT.*.bt2 $TMPDIR/

#path to reads fastq file
PATH_READS_DIRECTORY=#pathReadsDirectory
FASTQ_READ1=#read1
FASTQ_READ2=#read2

#output prefix
OUTPUT_PREFIX=#outputPrefix

#tophat output directory
PATH_OUTPUT_DIR=#pathOutputDir

#tophat parameters
MULT_READS=#multReads
EDIT_DIST0=#editDist0
             
echo "`${NOW}`copying reads to temporary scratch space..."
#check if mate file exists
if [ ! -f $PATH_READS_DIRECTORY/$FASTQ_READ2 ]
then

    echo "`{$NOW}`No mate file found. Skipped."   	
    exit 1

fi

#generate string that will replace empty and short sequence lines in fastq files
READ_LENGTH=`gzip -d -c $PATH_READS_DIRECTORY/$FASTQ_READ1 | head -n 100 | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq | tail -n 1`
STRING=$(for i in `eval echo {1..$READ_LENGTH}`;do printf "%s" "N";done;)

#replace empty sequence lines (together with corresponding quality lines) created by cutadapt
#unzip and copy fastqs to tmp space
FASTQ_READ1_NO_EXT=${FASTQ_READ1%.*}
echo "`${NOW}`$FASTQ_READ1_NO_EXT"
gzip -c -d $PATH_READS_DIRECTORY/$FASTQ_READ1 | sed "s/^$/$STRING/g" > $TMPDIR/$FASTQ_READ1_NO_EXT

FASTQ_READ2_NO_EXT=${FASTQ_READ2%.*}
echo "`${NOW}`$FASTQ_READ2_NO_EXT"
gzip -c -d $PATH_READS_DIRECTORY/$FASTQ_READ2 | sed "s/^$/$STRING/g" > $TMPDIR/$FASTQ_READ2_NO_EXT

COUNT_READS1=`wc -l $TMPDIR/$FASTQ_READ1_NO_EXT|cut -f 1 -d ' '`
COUNT_READS2=`wc -l $TMPDIR/$FASTQ_READ2_NO_EXT|cut -f 1 -d ' '`

#check if number of lines is the same in mate files
if [ "$COUNT_READS1" != "$COUNT_READS2" ]
then

    echo "`${NOW}`Unequal number of lines in the mate files. Skipped." 
    exit 1

fi

INSTRUMENT=`head -n 1 $TMPDIR/$FASTQ_READ1_NO_EXT|awk -F ':' '{print $1}'`
grep "$INSTRUMENT" $TMPDIR/$FASTQ_READ1_NO_EXT|cut -f 1 -d ' ' > $TMPDIR/R1.names
grep "$INSTRUMENT" $TMPDIR/$FASTQ_READ2_NO_EXT|cut -f 1 -d ' ' > $TMPDIR/R2.names
SYNC=`diff R1.names R2.names|head -n 1`

#check if mate files are sincronized
if [ ! -z $SYNC ]
then

    echo "`${NOW}`Mate files are not syncronized. Skipped." 
    exit 1

fi

#######################
#run TopHat alignment #
#######################

SECONDARY_ALIGNMENT_ARG=""
if [[ $MULT_READS == "T" ]]
then
	SECONDARY_ALIGNMENT_ARG="--report-secondary-alignments"
fi

EDIT_DIST0_ARG=""
if [[ $EDIT_DIST0 == "T" ]]
then
	EDIT_DIST0_ARG="--read-realign-edit-dist=0"
fi

echo "`${NOW}`running tophat alignment"
tophat $SECONDARY_ALIGNMENT_ARG $EDIT_DIST0_ARG -p $THREADS --transcriptome-index=$TMPDIR/$ANNOTATION_BASENAME $TMPDIR/$REFERENCE_FASTA_BASENAME $TMPDIR/$FASTQ_READ1_NO_EXT $TMPDIR/$FASTQ_READ2_NO_EXT

echo "`${NOW}`merging accepted_hits.bam and unmapped.bam"

java -jar -Xmx$JAVA_XMX $PICARD_HOME/MergeSamFiles.jar INPUT=$TMPDIR/tophat_out/accepted_hits.bam INPUT=$TMPDIR/tophat_out/unmapped.bam OUTPUT=$TMPDIR/tophat_out/all.bam SORT_ORDER=coordinate ASSUME_SORTED=false MERGE_SEQUENCE_DICTIONARIES=false USE_THREADING=true VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR VERBOSITY=WARNING
 
java -jar -Xmx$JAVA_XMX $PICARD_HOME/ReorderSam.jar INPUT=$TMPDIR/tophat_out/all.bam OUTPUT=$TMPDIR/tophat_out/all.sorted.bam REFERENCE=$TMPDIR/$REFERENCE_FASTA_BASENAME.fa VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR VERBOSITY=WARNING

# clean BAM (soft-clip an alignment that hangs off the end 
# of its reference sequence and set MAPQ to 0 if read is unmapped)
echo "`${NOW}`cleaning merged BAM..."
java -jar -Xmx$JAVA_XMX $PICARD_HOME/CleanSam.jar INPUT=$TMPDIR/tophat_out/all.sorted.bam OUTPUT=$TMPDIR/tophat_out/all.sorted.clean.bam VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR VERBOSITY=WARNING

# index output BAM
echo "`${NOW}`indexing merged BAM..."
samtools index $TMPDIR/tophat_out/all.sorted.clean.bam

echo "`${NOW}`copying results to $PATH_OUTPUT_DIR"
cp $TMPDIR/tophat_out/all.sorted.clean.bam $PATH_OUTPUT_DIR/$OUTPUT_PREFIX.sorted.bam
cp $TMPDIR/tophat_out/all.sorted.clean.bam.bai $PATH_OUTPUT_DIR/$OUTPUT_PREFIX.sorted.bam.bai
cp $TMPDIR/tophat_out/align_summary.txt $PATH_OUTPUT_DIR/$OUTPUT_PREFIX.align_summary.txt

mkdir $PATH_OUTPUT_DIR/$OUTPUT_PREFIX.logs
cp -r $TMPDIR/tophat_out/logs/* $PATH_OUTPUT_DIR/$OUTPUT_PREFIX.logs

chmod 660 $PATH_OUTPUT_DIR/*
chmod 770 $PATH_OUTPUT_DIR/*.logs
chmod 660 $PATH_OUTPUT_DIR/*.logs/*

ls -l 
ls -l $TMPDIR/tophat_out
ls -l $TMPDIR/tophat_out/logs

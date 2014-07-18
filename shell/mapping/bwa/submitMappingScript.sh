#!/bin/bash

#
# script to submit BWA mapping scripts after 
# fastq split is complete

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=800mb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

WALLTIME_HOURS_PER_RUN=#walltimeHoursPerRun

QUEUE=#queue
BASEDIR=#baseDir
SETUP_LOG=#setupLog
PATH_TMP_DIR=#pathTmpDir
FASTQ_READ1_NO_EXT=#fastqRead1NoExt
PATTERN_READ_1=#patternRead1
PATTERN_READ_2=#patternRead2
THREADS_PER_RUN=#threadsPerRun
PATH_SCRIPTS_DIR=#pathScriptsDir
PATH_RESULTS_DIR=#pathResultsDir
PATH_MAPPING_DIR=#pathMappingDir
PATH_REFERENCE_FASTA_NO_EXT=#pathReferenceFastaNoExt
REFERENCE_FASTA_NAME=`basename $PATH_REFERENCE_FASTA_NO_EXT .fa`
PATH_RUN_DIR=#pathRunDir
TODAY=#today
DEPLOYMENT_SERVER=#deploymentServer
SUMMARY_DEPLOYMENT=#summaryDeployment
SUMMARY_RESULTS=#summaryResults
TRUE_PROJECT_DIR=#isProjectDir

#variables to store job dependencies
MERGE_DEPENDENCIES=afterok

#variable to store files to merge
MERGE_FILES=""


#for each subset of reads...
echo "`$NOW`submitting mapping scripts:"
for FASTQ_READ1_SPLIT in `ls --color=never $PATH_TMP_DIR/${FASTQ_READ1_NO_EXT}_split/*.f*q* | grep -v $PATTERN_READ_2`
do

	#remove path information
	FASTQ_READ1_SPLIT=`basename $FASTQ_READ1_SPLIT`

	#get name of read2 fastq by preplacing read pair tag
	FASTQ_READ2_SPLIT=`echo $FASTQ_READ1_SPLIT | perl -pe "s/$PATTERN_READ_1/$PATTERN_READ_2/"`

	PATH_READS_FASTQ_READ1_SPLIT=$PATH_TMP_DIR/${FASTQ_READ1_NO_EXT}_split/$FASTQ_READ1_SPLIT
	PATH_READS_FASTQ_READ2_SPLIT=$PATH_TMP_DIR/${FASTQ_READ1_NO_EXT}_split/$FASTQ_READ2_SPLIT
     
	#output prefix
	OUTPUT_PREFIX=$FASTQ_READ1_SPLIT.vs.$REFERENCE_FASTA_NAME

	SCRIPT_PATH=$PATH_SCRIPTS_DIR/bwaAlignPe.$OUTPUT_PREFIX.sh
	cp $BASEDIR/bwaAlignPe.sh $SCRIPT_PATH
	chmod 770 $SCRIPT_PATH

	#set variables 
	sed -i -e "s/walltimeHours/$WALLTIME_HOURS_PER_RUN/" $SCRIPT_PATH  
	sed -i -e "s/threads/$THREADS_PER_RUN/" $SCRIPT_PATH
	sed -i -e "s/outputPrefix/$OUTPUT_PREFIX/" $SCRIPT_PATH
	sed -i -e "s/multReads/$MULT_READS/" $SCRIPT_PATH
	sed -i -e "s/pathOutputDir/${PATH_MAPPING_DIR//\//\\/}/" $SCRIPT_PATH
	sed -i -e "s/pathReferenceFastaNoExt/${PATH_REFERENCE_FASTA_NO_EXT//\//\\/}/" $SCRIPT_PATH
	sed -i -e "s/pathReadsFastqRead1NoExt/${PATH_READS_FASTQ_READ1_SPLIT//\//\\/}/" $SCRIPT_PATH
	sed -i -e "s/pathReadsFastqRead2NoExt/${PATH_READS_FASTQ_READ2_SPLIT//\//\\/}/" $SCRIPT_PATH

	#submit job and save job ID to dependency variable 
	LOG_OUTPUT_PATH=`echo $SCRIPT_PATH | perl -pe 's/\.sh/\.log/g'`
	echo "`$NOW`bwaAlignPe.$OUTPUT_PREFIX.sh"
	echo -n "`$NOW`"
	JOB_ID=`qsub -o $LOG_OUTPUT_PATH $SCRIPT_PATH`
	echo $JOB_ID
	MERGE_DEPENDENCIES=$MERGE_DEPENDENCIES:$JOB_ID 
	MERGE_FILES="$MERGE_FILES $OUTPUT_PREFIX.unsorted.bam"

done;
  
#submit merging jobs
   
OUTPUT_PREFIX=$FASTQ_READ1_NO_EXT.vs.$REFERENCE_FASTA_NAME

SCRIPT_PATH=$PATH_SCRIPTS_DIR/samtoolsMerge.$OUTPUT_PREFIX.sh
cp $BASEDIR/../samtools/samtoolsMergeAndDelete.sh $SCRIPT_PATH
chmod 770 $SCRIPT_PATH

sed -i -e "s/outputPrefix/$OUTPUT_PREFIX/" $SCRIPT_PATH
sed -i -e "s/inputDir/${PATH_MAPPING_DIR//\//\\/}/" $SCRIPT_PATH
sed -i -e "s/pathOutputDir/${PATH_RESULTS_DIR//\//\\/}/" $SCRIPT_PATH
sed -i -e "s/inBam/\"${MERGE_FILES//\//\\/}\"/" $SCRIPT_PATH

LOG_OUTPUT_PATH=`echo $SCRIPT_PATH | perl -pe 's/\.sh/\.log/g'`
echo "`$NOW`submitting merging script:"
echo "`$NOW`samtoolsMerge.$OUTPUT_PREFIX.sh"
echo -n "`$NOW`"

MERGE_JOB_ID=`qsub  -o $LOG_OUTPUT_PATH -W depend=$MERGE_DEPENDENCIES $SCRIPT_PATH`
echo $MERGE_JOB_ID

SUMMARY_SCRIPT=$PATH_SCRIPTS_DIR/summary_bwa.$OUTPUT_PREFIX.pl
cp $BASEDIR/summary_bwa.pl $SUMMARY_SCRIPT
chmod 770 $SUMMARY_SCRIPT

sed -i -e "s/isProject/$TRUE_PROJECT_DIR/" $SUMMARY_SCRIPT
sed -i -e "s/projectDirAnalysis/${PATH_RUN_DIR//\//\\/}/" $SUMMARY_SCRIPT
sed -i -e "s/projectDirResults/${PATH_RESULTS_DIR//\//\\/}/" $SUMMARY_SCRIPT
sed -i -e "s/Today/$TODAY/" $SUMMARY_SCRIPT
sed -i -e "s/deploymentServer/$DEPLOYMENT_SERVER/" $SUMMARY_SCRIPT
sed -i -e "s/summaryDeployment/${SUMMARY_DEPLOYMENT//\//\\/}/" $SUMMARY_SCRIPT
sed -i -e "s/summaryResults/${SUMMARY_RESULTS//\//\\/}/" $SUMMARY_SCRIPT

SUM_DEPENDENCIES=afterany:$MERGE_JOB_ID
SUMMARY_LOG=`echo $SUMMARY_SCRIPT | perl -pe 's/\.pl/\.log/g'`
echo "`$NOW`submitting summary script:"
echo "`$NOW`$SUMMARY_SCRIPT"
echo -n "`$NOW`"

SUM_JOB_ID=`qsub -q $QUEUE -o $SUMMARY_LOG  -j oe -W depend=$SUM_DEPENDENCIES -M cgi@imperial.ac.uk $SUMMARY_SCRIPT`
echo $SUM_JOB_ID


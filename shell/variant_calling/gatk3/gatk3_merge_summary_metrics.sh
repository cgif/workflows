#!/bin/bash

## script to merge and plot GATK coverage metrics 

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=8gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules
module load R/#rVersion

RESULTS_DIR=#resultsDir
RUN_DIR=#runDir
TYPE=#type
TRANSPOSE_TABLE_SCRIPT=#transposeTableScript

cp $TRANSPOSE_TABLE_SCRIPT $TMPDIR/transposeTable.R

#RESULTS_DIR=$1
#RUN_DIR=$2
#TYPE=$3

NOW="date +%Y-%m-%d%t%T%t"

SCRIPT_CODE="GATKMSME"

RUN_LOG=#runLog

LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"

RESULTS_DIR_RUN=`dirname $RESULTS_DIR`

RUN_DATE=`basename $RESULTS_DIR_RUN`

BASEDIR="$( cd "$( dirname "$0" )" && pwd )"

#get project name
PROJECT_NAME=`echo $RESULTS_DIR_RUN | perl -e 'while(<>){ @cols=split(/\//); print $cols[4]; }'`


#make output directory
#mkdir -p $RESULTS_DIR/metrics

#merge cumulative coverage and sample summaries

for METRIC in sample_cumulative_coverage_proportions sample_summary
do


	echo "`${NOW}`INFO $SCRIPT_CODE merging $METRIC metrics..."
	#get first file to extract header from
	HEADER_FILE=`ls $RESULTS_DIR_RUN/*/recalibration/*.bam.coverage.$METRIC | head -n 1`
	OUTPUT_FILE=$RESULTS_DIR_RUN/multisample/metrics/$PROJECT_NAME.$RUN_DATE.$METRIC
	OUTPUT_FILE_TMP=$TMPDIR/$METRIC
	
	head -n 1 $HEADER_FILE > $OUTPUT_FILE_TMP
	
	for SAMPLE in `ls $RESULTS_DIR_RUN`
	do

		if [[ "$SAMPLE" != "multisample" ]]
		then

			#echo $SAMPLE
			METRIC_FILE=$RESULTS_DIR_RUN/$SAMPLE/recalibration/$SAMPLE.bam.coverage.$METRIC
			if [[ -f "$METRIC_FILE" ]]
			then
				head -n 2 $METRIC_FILE | tail -n 1  >> $OUTPUT_FILE_TMP
			else
				#echo $METRIC_FILE
				echo "Warning: $METRIC_FILE does not exist"
			fi

		fi
		
	done	

	cp $OUTPUT_FILE_TMP $OUTPUT_FILE

	#logging
	STATUS=OK
	if [[ ! -s $OUTPUT_FILE ]]
	then
		STATUS=FAILED
	fi

	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t${METRIC}_metrics_merged\t$STATUS" >> $RUN_LOG

	echo "`${NOW}`INFO $SCRIPT_CODE done..."

done

#merge interval summary
METRIC=sample_interval_summary

HEADER_FILE=`ls $RESULTS_DIR_RUN/*/recalibration/*.bam.coverage.$METRIC | head -n 1`
OUTPUT_FILE=$RESULTS_DIR_RUN/multisample/metrics/$PROJECT_NAME.$RUN_DATE.$METRIC
OUTPUT_FILE_TMP=$TMPDIR/$METRIC

#print header
cut -f1 $HEADER_FILE | perl -e 'while(<>){ s/Target/Sample/; if(!eof){s/\n/\t/}; print;  }' > $OUTPUT_FILE_TMP

echo "`${NOW}`INFO $SCRIPT_CODE merging $METRIC metrics..."

for SAMPLE in `ls $RESULTS_DIR_RUN`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then

		#echo $SAMPLE
		METRIC_FILE=$RESULTS_DIR_RUN/$SAMPLE/recalibration/$SAMPLE.bam.coverage.$METRIC
		if [[ -f "$METRIC_FILE" ]]
		then
			cut -f3 $METRIC_FILE | perl -e "while(<>){ 
												if(/average_coverage/){
													print \"$SAMPLE\t\"; 
												} else {
													s/\n//; 
													printf(\"%.0f\", \$_); 
													if(!eof){ 
														print \"\t\"; 
													} else { 
														print \"\n\"; 
													} 
												}
											}" >> $OUTPUT_FILE_TMP
		else
			echo "Warning: $METRIC_FILE does not exist"
		fi

	fi
		
done	

echo "`${NOW}`INFO $SCRIPT_CODE done"

#transpose interval summary

echo "`${NOW}`INFO $SCRIPT_CODE transposing interval_summary table..."

#configure script
sed -i -e "s/#inputTable/${OUTPUT_FILE_TMP//\//\\/}/" $TMPDIR/transposeTable.R

#execute script
R --vanilla < $TMPDIR/transposeTable.R

#reformat chromosome coordinates
perl -i -pe 's/X(.?)\.(\d*?)\.(\d*?)/$1:$2-$3/' $OUTPUT_FILE_TMP.tr

mv $OUTPUT_FILE_TMP.tr $OUTPUT_FILE_TMP

cp $OUTPUT_FILE_TMP $OUTPUT_FILE

#logging
STATUS=OK
if [[ ! -s $OUTPUT_FILE ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t${METRIC}_metrics_merged\t$STATUS" >> $RUN_LOG


echo "`${NOW}`INFO $SCRIPT_CODE done"


#generate XLSX from TSV

echo "`${NOW}`INFO $SCRIPT_CODE generating Excel spreadsheets..."

for METRIC in sample_cumulative_coverage_proportions sample_summary sample_interval_summary
do

	INPUT_FILE=$TMPDIR/$METRIC
	OUTPUT_FILE=$RESULTS_DIR_RUN/multisample/metrics/$PROJECT_NAME.$RUN_DATE.$METRIC.xlsx

	$RUN_DIR/tsvToXls.pl $INPUT_FILE $OUTPUT_FILE

	#logging
	STATUS=OK
	if [[ ! -s $OUTPUT_FILE ]]
	then
		STATUS=FAILED
	fi

	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t${METRIC}_metrics_merged_xlsx\t$STATUS" >> $RUN_LOG


done

echo "`${NOW}`INFO $SCRIPT_CODE done"


#plot coverage across samples and amplicons
echo "`${NOW}`INFO $SCRIPT_CODE generating coverage plots..."

if [[ "$TYPE" == "TARGETED" ]]; then
	
	R --vanilla < $RUN_DIR/plot_summary_metrics.R
	
fi
echo "`${NOW}`INFO $SCRIPT_CODE done"


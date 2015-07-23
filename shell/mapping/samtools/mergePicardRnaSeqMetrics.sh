#!/bin/bash

## script to merge RnaSeqMetrics files

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=800mb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

MERGETAG_PROJECT_DIRECTORY=mergeTagProjectDirectory
MERGETAG_DATE=mergeTagDate
PROJECT_NAME=mergeTagProjectName

OUTPUT_DIR=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/multisample

echo "`$NOW`merging RnaSeqMetrics files..."

#get first file to extract header from
HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.RnaSeqMetrics | head -n 1`

#create outputfile
OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.RnaSeqMetrics

#extract column headers and append them to output file
echo -n -e "#SAMPLE\t" > $OUTPUT_FILE
grep 'PF_BASES' $HEADER_FILE | perl -pe 's/SAMPLE\tLIBRARY\tREAD_GROUP\n/\n/' >> $OUTPUT_FILE

#create outputfile for RNA integrity charts
OUTPUT_CHART=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.chartOutput.pdf
HEADER_CHART=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.chartOutput | head -n 1`
cp $HEADER_CHART $OUTPUT_CHART

#append sample RnaSeq metrics
for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then
	
		METRICS_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.RnaSeqMetrics
		if [[ -e "$METRICS_FILE" ]]		
		then		

			echo "`$NOW`$SAMPLE"
			#extract stats, remove blank lines and fill in sample information

			cat $METRICS_FILE | perl -e "while(<>){print \"$SAMPLE\t\$_\" if /^\d+/; last if /^\d+/;} " >> $OUTPUT_FILE
		
		fi

		CHART_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.chartOutput
		if [[ -e "$CHART_FILE" ]]		
		then		

		        if [[ "$CHART_FILE" != "$HEADER_CHART" ]]
			then

			        #merge pdf images for RNA integrity charts
			        convert $OUTPUT_CHART $CHART_FILE $OUTPUT_CHART

			fi
		
		fi
		
	fi

done

#plot RS metrics
echo "`${NOW}`generating plots for RS metrics..."

R --vanilla < $RUN_DIR/plot_RS_metrics.R


#!/bin/bash

## script to merge AlignmentSummaryMetrics files

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=800mb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

MERGETAG_PROJECT_DIRECTORY=mergeTagProjectDirectory
MERGETAG_DATE=mergeTagDate
PROJECT_NAME=mergeTagProjectName

OUTPUT_DIR=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/multisample

echo "`$NOW`merging AlignmentSummaryMetrics files..."

#get first file to extract header from
HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*alignment_summary_metrics | head -n 1`

for CATEGORY in FIRST_OF_PAIR SECOND_OF_PAIR PAIR UNPAIRED
do

	#create outputfile
	OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.alignment_summary_metrics.$CATEGORY

	#extract column headers and append them to output file
	echo -n "#" > $OUTPUT_FILE
	grep 'CATEGORY' $HEADER_FILE | perl -pe 's/CATEGORY/SAMPLE/' | perl -pe 's/SAMPLE\tLIBRARY\tREAD_GROUP\n/\n/' >> $OUTPUT_FILE

done;

#append sample alignment summary metrics
for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then
	
		METRICS_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.alignment_summary_metrics
		if [[ -e "$METRICS_FILE" ]]		
		then		
	
			for CATEGORY in FIRST_OF_PAIR SECOND_OF_PAIR PAIR UNPAIRED
			do

				OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.alignment_summary_metrics.$CATEGORY
				
				#extract stats, remove blank lines and fill in amplicon set, sample, library and readgroup information
				grep -P "^$CATEGORY\t" $METRICS_FILE | \
					perl -e "while(<>){ chomp;
								@cols=split(/\t/); 
								print \"$SAMPLE\";
								for (\$i=1; \$i < @cols; \$i++){
									
									print \"\t\$cols[\$i]\";
	
								}
								print \"\n\";
							} " \
					>> $OUTPUT_FILE
		
			done

		fi
		
	fi

done


#!/bin/bash

## script to merge HybridSequencingMetrics files

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
CUSTOM_AMPLICON_SET=customAmpliconSet
RUN_DIR=#runFolder

OUTPUT_DIR=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/multisample

echo "`$NOW`merging hybridMetrics files..."

#get first file to extract header from
HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.bam.hybridMetrics | head -n 1`

#create outputfile
OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.hybridMetrics

#extract column headers and append them to output file
echo -n -e "#SAMPLE\t" > $OUTPUT_FILE
grep 'BAIT_SET' $HEADER_FILE | perl -pe 's/SAMPLE\tLIBRARY\tREAD_GROUP\n/\n/' >> $OUTPUT_FILE

#append sample targeted PCR metrics
for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then
	
		METRICS_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.bam.hybridMetrics
		if [[ -e "$METRICS_FILE" ]]		
		then		
			echo "`$NOW`$SAMPLE"
			#extract stats, remove blank lines and fill in amplicon set, sample, library and readgroup information
			grep -vP '#|BAIT_SET' $METRICS_FILE | grep -v '^$' | \
				perl -e "while(<>){ chomp;
							@cols=split(/\t/); 
							print \"$SAMPLE\t$CUSTOM_AMPLICON_SET\";
							for (\$i=1; \$i < @cols; \$i++){
								
								print \"\t\$cols[\$i]\";

							}
							print \"\n\";
						} " \
				>> $OUTPUT_FILE
		
		fi
		
	fi

done


echo "`$NOW`merging perTargetCoverage mean_coverage information..."
#merge perTargetCoverage
HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.bam.perTargetCoverage | head -n 1`

OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.perTargetCoverage

#append column headers with header
echo -n -e "#SAMPLE" > $OUTPUT_FILE
sed 1d $HEADER_FILE | perl -e 'while(<>){ @cols=split(/\t/); print "\t".$cols[4]."[".$cols[0].":".$cols[1]."-".$cols[2]."]"; } print "\n";' >> $OUTPUT_FILE

for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then
		
		METRICS_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.bam.perTargetCoverage
		if [[ -e "$METRICS_FILE" ]]		
		then
				
			echo "`$NOW`$SAMPLE"
			echo -n $SAMPLE >> $OUTPUT_FILE
			sed 1d $METRICS_FILE | perl -e 'while(<>){ @cols=split(/\t/); print "\t"; printf("%.0f", $cols[6]); } print "\n";' >> $OUTPUT_FILE

		fi
		
	fi
	
done

#plot HS metrics
echo "`${NOW}`generating plots for HS metrics..."

R --vanilla < $RUN_DIR/plot_HS_metrics.R

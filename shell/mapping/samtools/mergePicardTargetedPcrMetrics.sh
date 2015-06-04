#!/bin/bash

## script to merge targeted PCR metrics

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=800mb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

module load R/3.1.0

NOW="date +%Y-%m-%d%t%T%t"

MERGETAG_PROJECT_DIRECTORY=mergeTagProjectDirectory
MERGETAG_DATE=mergeTagDate
PROJECT_NAME=mergeTagProjectName
CUSTOM_AMPLICON_SET=customAmpliconSet
NON_OVERLAPPING=#nonOverlapping
RUN_DIR=#runFolder


OUTPUT_DIR=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/multisample

echo "`$NOW`merging targetedPcrMetrics files..."

#get first file to extract header from
HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.bam.targetedPcrMetrics | head -n 1`

#create outputfile
OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.targetedPcrMetrics

#extract column headers and append them to output file
echo -n -e "#SAMPLE\t" > $OUTPUT_FILE
grep 'CUSTOM_AMPLICON_SET' $HEADER_FILE | perl -pe 's/SAMPLE\tLIBRARY\tREAD_GROUP\n/\n/' >> $OUTPUT_FILE

#append sample targeted PCR metrics
for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then
	
		METRICS_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.bam.targetedPcrMetrics
		if [[ -e "$METRICS_FILE" ]]		
		then		
			echo "`$NOW`$SAMPLE"
			#extract stats, remove blank lines and fill in amplicon set, sample, library and readgroup information
			grep -vP '#|CUSTOM_AMPLICON_SET' $METRICS_FILE | grep -v '^$' | \
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

#extract column headers and append them to output file
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

# if we have non-overlapping intervals statistics, merge it as well.

if [[ "$NON_OVERLAPPING" == TRUE ]]; then
	
	echo "`$NOW`merging perTargetCoverage mean_coverage information for non-overlapping amplicon regions..."
	#merge non_overlapping.perTargetCoverage
	HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.bam.non_overlapping.perTargetCoverage | head -n 1`

	OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.non_overlapping.perTargetCoverage

	#extract column headers and append them to output file
	echo -n -e "#SAMPLE" > $OUTPUT_FILE
	sed 1d $HEADER_FILE | perl -e 'while(<>){ @cols=split(/\t/); print "\t".$cols[4]."[".$cols[0].":".$cols[1]."-".$cols[2]."]"; } print "\n";' >> $OUTPUT_FILE

	for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`; do

		if [[ "$SAMPLE" != "multisample" ]]; then
		
			METRICS_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.bam.non_overlapping.perTargetCoverage
			if [[ -e "$METRICS_FILE" ]]; then
				
				echo "`$NOW`$SAMPLE"
				echo -n $SAMPLE >> $OUTPUT_FILE
				sed 1d $METRICS_FILE | perl -e 'while(<>){ @cols=split(/\t/); print "\t"; printf("%.0f", $cols[6]); } print "\n";' >> $OUTPUT_FILE

			fi
		fi
	done
fi

#merge interval coverage from GATK output 

echo "`${NOW}` merging interval coverage from GATK DepthOfCoverage..."
METRIC=sample_interval_summary

HEADER_FILE=`ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/*/*.bam.perAmpliconCoverage.$METRIC | head -n 1`
OUTPUT_FILE=$OUTPUT_DIR/$PROJECT_NAME.$MERGETAG_DATE.$METRIC

#print header
cut -f1 $HEADER_FILE | perl -e 'while(<>){ s/Target/Sample/; if(!eof){s/\n/\t/}; print;  }' > $OUTPUT_FILE

for SAMPLE in `ls $MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE`
do

	if [[ "$SAMPLE" != "multisample" ]]
	then

		METRIC_FILE=$MERGETAG_PROJECT_DIRECTORY/$MERGETAG_DATE/$SAMPLE/$SAMPLE.bam.perAmpliconCoverage.$METRIC
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
											}" >> $OUTPUT_FILE
		else
			echo "Warning: $METRIC_FILE for $SAMPLE does not exist"
		fi

	fi
		
done	

perl -i -pe 's/X(.*?)\.(\d*?)\.(\d*?)/$1:$2-$3/' $OUTPUT_FILE

#plot coverage across samples and amplicons
echo "`${NOW}`generating coverage plots..."

R --vanilla < $RUN_DIR/plot_amplicon_summary_metrics.R
	
echo "`${NOW}`done"



#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=800mb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load java/#javaVersion

NXTGENUTILS_VERSION=#nxtGenUtilsVersion
NXTGENUTILS_HOME=/project/tgu/bin/nxtgen-utils-${NXTGENUTILS_VERSION}


RECALIBRATION_REPORTS="#recalibrationReports"
MERGED_RECALIBRATION_REPORT=#mergedRecalibrationReport
RUN_LOG=#runLog
SAMPLE=#sample

SCRIPT_CODE="GATKMERR"

LOG_INFO="`${NOW}`INFO $SCRIPT_CODE"
LOG_ERR="`${NOW}`ERROR $SCRIPT_CODE"
LOG_WARN="`${NOW}`WARN $SCRIPT_CODE"
LOG_DEBUG="`${NOW}`DEBUG $SCRIPT_CODE"


echo "`${NOW}`INFO $SCRIPT_CODE copying recalibration reports to tmp directory..."
INPUT_REPORT=""
REPORT_COUNT=0
for REPORT in $RECALIBRATION_REPORTS; do
	
	REPORT_COUNT=$(( $REPORT_COUNT + 1 ))
	TMP_REPORT=$TMPDIR/in_report_$REPORT_COUNT	
	cp $REPORT $TMP_REPORT
    INPUT_REPORT="$INPUT_REPORT -i $TMP_REPORT"

done

# step 5: merge recalibration reports
echo "`${NOW}`INFO $SCRIPT_CODE merging recalibration reports..."
echo "`${NOW}`DEBUG $SCRIPT_CODE java -jar nxtgen-utils-0.12.3/NxtGenUtils.jar GatherGatkBqsrReports $INPUT_REPORT -o $TMPDIR/merged_recal_data.grp"
java -jar $NXTGENUTILS_HOME/NxtGenUtils.jar GatherGatkBqsrReports $INPUT_REPORT -o $TMPDIR/merged_recal_data.grp

echo "`${NOW}`INFO $SCRIPT_CODE copying merged recalibration report to $MERGED_RECALIBRATION_REPORT..."
cp $TMPDIR/merged_recal_data.grp $MERGED_RECALIBRATION_REPORT


#logging
STATUS=OK
if [[ ! -e $MERGED_RECALIBRATION_REPORT ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tprerecalibration_report_merged\t$STATUS" >> $RUN_LOG


echo "`${NOW}`INFO $SCRIPT_CODE done"


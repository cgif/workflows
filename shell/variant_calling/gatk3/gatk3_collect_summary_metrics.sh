#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=10gb:tmpspace=#tmpSpaceMbmb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load gatk/#gatkVersion
module load R/#rVersion
module load java/#javaVersion

#NxtGenUtils path
NXTGENUTILS_VERSION=#nxtGenUtilsVersion
NXTGENUTILS_HOME=/groupvol/cgi/bin/nxtgen-utils-${NXTGENUTILS_VERSION}

#Java max. memory
JAVA_XMX=4G

SCRIPT_CODE="GATKCSME"

LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"


#REALIGNED_RECALIBRATED_BAM=#reaglignedRecalibratedBam
SAMPLE=#sample
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
ANALYSIS_DIR=#analysisDir
RESULTS_DIR=#resultsDir
REALIGNED_RECALIBRATED_BAM=$RESULTS_DIR/recalibration/$SAMPLE.bam
TYPE=#sequencingType

MERGED_PRE_RECALIBRATION_REPORT=#mergedPreRecalibrationReport
POST_RECALIBRATION_REPORTS="#recalibrationReports"
MERGED_POST_RECALIBRATION_REPORT=#mergedPostRecalibrationReport
POST_RECALIBRATION_PLOTS_OUTPUT_DIR=#postRecalibrationPlotsOutputDir

TARGET_INTERVALS_FILE=#targetIntervals
SUMMARY_SCRIPT_PATH=#summaryScriptPath
RUN_LOG=#runLog


echo "`${NOW}`copying reference fasta and indexto tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

echo "`${NOW}`copying realigned BAM to tmp directory..."
BAM_NAME=`basename $REALIGNED_RECALIBRATED_BAM`
BAM_DIR=`dirname $REALIGNED_RECALIBRATED_BAM`
cp $REALIGNED_RECALIBRATED_BAM $TMPDIR/$BAM_NAME
cp $REALIGNED_RECALIBRATED_BAM.bai $TMPDIR/$BAM_NAME.bai

echo "`${NOW}`copying merged pre-recalibration report to tmp directory..."
cp $MERGED_PRE_RECALIBRATION_REPORT $TMPDIR/merged_pre_recal_data.grp

echo "`${NOW}`copying chunk post-recalibration reports to tmp directory..."
INPUT_REPORT=""
REPORT_COUNT=0
for REPORT in $POST_RECALIBRATION_REPORTS; do
	
	REPORT_COUNT=$(( $REPORT_COUNT + 1 ))
	TMP_REPORT=$TMPDIR/in_report_$REPORT_COUNT	
	cp $REPORT $TMP_REPORT
        INPUT_REPORT="$INPUT_REPORT -i $TMP_REPORT"

done

#merge post-recalibration reports
echo "`${NOW}`merging recalibration reports..."
echo "`${NOW}`java -jar nxtgen-utils-0.12.3/NxtGenUtils.jar GatherGatkBqsrReports $INPUT_REPORT -o $TMPDIR/merged_recal_data.grp"
java -jar $NXTGENUTILS_HOME/NxtGenUtils.jar GatherGatkBqsrReports $INPUT_REPORT -o $TMPDIR/merged_post_recal_data.grp

echo "`${NOW}`copying merged recalibration report to $MERGED_RECALIBRATION_REPORT..."
cp $TMPDIR/merged_post_recal_data.grp $MERGED_POST_RECALIBRATION_REPORT


#logging
STATUS=OK
if [[ ! -s $MERGED_POST_RECALIBRATION_REPORT ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\tpostrecalibration_report_merged\t$STATUS" >> $RUN_LOG

#generate post-recalibration plots

INTERVAL_ARG=""
if [ "$TARGET_INTERVALS_FILE" != "NULL" ]
then
	INTERVAL_ARG="-L $TARGET_INTERVALS_FILE"
fi 

echo "`${NOW}`generating recalibration plots for realigned and recalibrated BAM..."

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
    -T AnalyzeCovariates \
    -R $TMPDIR/reference.fa \
    -before $TMPDIR/merged_pre_recal_data.grp \
    -after $TMPDIR/merged_post_recal_data.grp \
    -plots $TMPDIR/${SAMPLE}.realigned.recalibrated.recal_plot.pdf \
    $INTERVAL_ARG

cp $TMPDIR/${SAMPLE}.realigned.recalibrated.recal_plot.pdf $POST_RECALIBRATION_PLOTS_OUTPUT_DIR/
convert $TMPDIR/${SAMPLE}.realigned.recalibrated.recal_plot.pdf $TMPDIR/${SAMPLE}.realigned.recalibrated.recal_plot.jpeg
cp $TMPDIR/${SAMPLE}.realigned.recalibrated.recal_plot-1.jpeg $POST_RECALIBRATION_PLOTS_OUTPUT_DIR/

#logging
STATUS=OK
if [[ ! -s $POST_RECALIBRATION_PLOTS_OUTPUT_DIR/${SAMPLE}.realigned.recalibrated.recal_plot.pdf ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\trecalibration_plot\t$STATUS" >> $RUN_LOG


#calculate depth of coverage
echo "`${NOW}`calculating depth of coverage"

#set coverage statistics accumulation
#categories
COVERAGE_ARG="-ct 2 -ct 4 -ct 10 -ct 20 -ct 30" 
if [[ "$TYPE" == "TARGETED" ]]; then
	COVERAGE_ARG="-ct 10 -ct 50 -ct 100 -ct 500 -ct 1000" 
fi


java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR -jar $GATK_HOME/GenomeAnalysisTK.jar \
  -T DepthOfCoverage \
  -R $REFERENCE_FASTA \
  -I $TMPDIR/$BAM_NAME \
  -o $TMPDIR/${BAM_NAME}.coverage \
  $COVERAGE_ARG \
  -omitBaseOutput \
  --countType COUNT_FRAGMENTS \
  -mbq 10 \
  -rf BadCigar \
  $INTERVAL_ARG

cp ${BAM_NAME}.coverage* $BAM_DIR/

#delete temporary files in analysis folder
#commented out because not tested properly yet
#rm $ANALYSIS_DIR/chunks/*bam
#rm $ANALYSIS_DIR/chunks/*bai
#rm $ANALYSIS_DIR/realignment/*bam
#rm $ANALYSIS_DIR/realignment/*bai
#rm $ANALYSIS_DIR/recalibration/*bam
#rm $ANALYSIS_DIR/recalibration/*bai

find $ANALYSIS_DIR -type f | xargs chmod 750 
find $RESULTS_DIR -type f | xargs chmod 750



## would be good to add -geneList refSeq.sorted.txt , but need to generate this list first 
echo "`${NOW}`done"


#diagnose targets
#TODO
#java -jar GenomeAnalysisTK.jar
#        -T DiagnoseTargets \
#              -R reference.fasta \
#              -o output.vcf \
#              -I sample1.bam \
#              -I sample2.bam \
#              -I sample3.bam \
#              -L intervals.interval_list


#run summary script
perl $SUMMARY_SCRIPT_PATH

#ls -alh
#ls -alh tmp


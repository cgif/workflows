#!/bin/bash

## script to run GATK for counting covariates before base quality recalibration

#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=#cpuThreads:mem=25gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion

NXTGENUTILS_HOME=/groupvol/cgi/bin/nxtgen-utils-#nextGenUtilsVersion

JAVA_XMX=24G
NCT=#nctThreads
RUN_LOG=#runLog
WARNING_LOG=#warningLog

SCRIPT_CODE="GATKPRRE"


LOG_INFO="`$NOW`INFO $SCRIPT_CODE"
LOG_ERR="`$NOW`ERROR $SCRIPT_CODE"
LOG_WARN="`$NOW`WARN $SCRIPT_CODE"
LOG_DEBUG="`$NOW`DEBUG $SCRIPT_CODE"


# define variables
REFERENCE_FASTA=#referenceFasta
SAMPLE=#sample
ANALYSIS_DIR=#analysisDir
REALIGNED_BAM_FILE=#realignedBamFile
RECALIBRATION_REPORT=#recalibrationReport
RECALIBRATED_OUTPUT_BAM=#recalibratedBamOutput
FRAGMENT_FILE=#fragmentFile
FRAGMENT=#fragmentName
INCLUDES_UNMAPPED=#includesUnmapped
SUMMARY_SCRIPT_PATH=#summaryScriptPath
TYPE=#type
CLIP_CYCLES=#clipCycles


mkdir $TMPDIR/tmp

echo "`${NOW}`INFO $SCRIPT_CODE copying realigned BAM and index files to tmp directory..."
cp $REALIGNED_BAM_FILE $TMPDIR/tmp.bam
cp $REALIGNED_BAM_FILE.bai $TMPDIR/tmp.bam.bai 

echo "`${NOW}`INFO $SCRIPT_CODE copying merged recalibration report to tmp directory..."
cp $RECALIBRATION_REPORT $TMPDIR/merged_recal_data.grp

echo "`${NOW}`INFO $SCRIPT_CODE copying reference to tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict

cp $FRAGMENT_FILE $TMPDIR/fragment.intervals

INTERVAL_ARG="-L $TMPDIR/fragment.intervals"
if [ "$INCLUDES_UNMAPPED" == "T" ]
then
	INTERVAL_ARG_PR="$INTERVAL_ARG -L unmapped"
else
	INTERVAL_ARG_PR=$INTERVAL_ARG
fi

#############################################
# print reads with recalibrated base calls
############################################
echo "`${NOW}`INFO $SCRIPT_CODE recalibrating chunk BAM..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
   -nct $NCT \
   -T PrintReads \
   -R $TMPDIR/reference.fa \
   -I $TMPDIR/tmp.bam \
   -BQSR $TMPDIR/merged_recal_data.grp \
   -o $TMPDIR/realigned.recalibrated.bam \
   $INTERVAL_ARG_PR

# samtools index recalibrated BAM file
echo "`${NOW}`INFO $SCRIPT_CODE indexing recalibrated BAM..."
samtools index $TMPDIR/realigned.recalibrated.bam

# clip primer sequences from amplicon sequencing
# DO NOT use option -CR WRITE_NS_Q0S, because it is messes up all quality scores
# GATK will ignor N bases, so it is enough to recode into Ns
# http://gatkforums.broadinstitute.org/discussion/3269/in-gatk-haplotype-caller-how-do-you-treat-n-the-undetermined-base-which-could-be-any-of-acgt
if [[ "$TYPE" == "TARGETED" ]] && [[ "$CLIP_CYCLES" -gt 0 ]]; then

	echo "`${NOW}`INFO $SCRIPT_CODE clipping first $CLIP_CYCLES from amplicon reads..."
	java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-T ClipReads \
	-R $TMPDIR/reference.fa \
	-I $TMPDIR/realigned.recalibrated.bam \
	-o $TMPDIR/realigned.recalibrated.clipped.bam \
	-CT "1-$CLIP_CYCLES" \
	-CR WRITE_NS \
	$INTERVAL_ARG_PR

	mv $TMPDIR/realigned.recalibrated.clipped.bam $TMPDIR/realigned.recalibrated.bam
	echo "`${NOW}`INFO $SCRIPT_CODE indexing recalibrated clipped BAM..."
	samtools index $TMPDIR/realigned.recalibrated.bam

fi

#check that realigned and recalibrated bams have the same number of reads, 
# if so, copy to destination folder

READ_COUNT_INPUT=`samtools flagstat $TMPDIR/tmp.bam | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
READ_COUNT_OUTPUT=`samtools flagstat $TMPDIR/realigned.recalibrated.bam | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`
echo "Number of reads in realigned file: $READ_COUNT_INPUT"
echo "Number of reads in recalibrated file: $READ_COUNT_OUTPUT"

echo "`${NOW}`INFO $SCRIPT_CODE copying recalibrated BAM to output directory $ANALYSIS_DIR..."
cp $TMPDIR/realigned.recalibrated.bam $RECALIBRATED_OUTPUT_BAM
cp $TMPDIR/realigned.recalibrated.bam.bai $RECALIBRATED_OUTPUT_BAM.bai


if [[ $READ_COUNT_INPUT -ne $READ_COUNT_OUTPUT ]]; then
	echo "WARNING: $SCRIPT_CODE realined and recalibrated bams have different numbers of reads"	
	echo "$SAMPLE $FRAGMENT $SCRIPT_CODE input realined bam and recalibrated bam have different numbers of reads" >> $WARNING_LOG
	echo "$SAMPLE $FRAGMENT $SCRIPT_CODE Number of reads in realined chunk $FRAGMENT bam file: $READ_COUNT_INPUT" >> $WARNING_LOG
	echo "$SAMPLE $FRAGMENT $SCRIPT_CODE Number of reads in recalibrated chunk $FRAGMENT bam file: $READ_COUNT_OUTPUT" >> $WARNING_LOG
fi

#logging and deleting intermediate file

if [[ -s $RECALIBRATED_OUTPUT_BAM ]] && [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]] ; then
	STATUS=OK
	echo "deleting realigned bam file $REALIGNED_BAM_FILE"
	rm $REALIGNED_BAM_FILE $REALIGNED_BAM_FILE.bai
elif [[ -s $RECALIBRATED_OUTPUT_BAM ]] && [[ $READ_COUNT_INPUT -ne $READ_COUNT_OUTPUT ]]; then
	STATUS=WARNING
	echo "WARNING: $SCRIPT_CODE keeping input file for checks $REALIGNED_BAM_FILE"
else
	STATUS=FAILED
	echo "ERROR: $SCRIPT_CODE recalibration failed, keeping realigned file for rerun $REALIGNED_BAM_FILE"
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\trecalibrated_bam\t$STATUS" >> $RUN_LOG

if  [[ "$STATUS" == "FAILED" ]]; then 
	echo "`${NOW}` realignement failed, exiting"
	exit 1
fi

#run summary script
perl $SUMMARY_SCRIPT_PATH

echo "`${NOW}`INFO $SCRIPT_CODE done"


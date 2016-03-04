#!/bin/bash

## script to run GATK

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=#dataThreads:mem=5gb:tmpspace=#tmpSpaceMbmb

#PBS -M igf@imperial.ac.uk
#PBS -m ae
#PBS -j oe

#PBS -q pqcgi

# load modules
module load gatk/#gatkVersion
module load samtools/#samtoolsVersion
#Picard Version needs to be > 1.85 for ClipPrimerSequences to
#work properly 
#module load picard/#picardVersion
module load R/#rVersion
module load java/#javaVersion

NXTGENUTILS_VERSION=#nxtGenUtilsVersion
NXTGENUTILS_HOME=/project/tgu/bin/nxtgen-utils-${NXTGENUTILS_VERSION}

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4800M

SCRIPT_CODE="GATKRARC"

RUN_LOG=#runLog
WARNING_LOG=#warningLog

LOG_INFO="`${NOW}`INFO $SCRIPT_CODE"
LOG_ERR="`${NOW}`ERR $SCRIPT_CODE"
LOG_WARN="`${NOW}`WARN $SCRIPT_CODE"
LOG_DEBUG="`${NOW}`DEBUG $SCRIPT_CODE"

# define variables

INPUT_BAM=#inputBam
INPUT_BAM_NAME=`basename $INPUT_BAM .bam`
REFERENCE_FASTA=#referenceFasta
REFRENCE_SEQ_DICT=`echo $REFERENCE_FASTA | perl -pe 's/\.fa/\.dict/'`
ANALYSIS_DIR=#analysisDir
FRAGMENT=#fragmentName
RTC_DATA_THREADS=#rtcDataThreads
SAMPLE=#sample
INCLUDES_UNMAPPED=#includesUnmapped
FRAGMENT_FILE=#fragmentFile
TARGET_FILE=#targetFile
PRIMER_COORD_BED=#primerCoordBed
PRIMER_COORD_OFFSET=#primerCoordOffset
SUMMARY_SCRIPT_PATH=#summaryScriptPath

cp $FRAGMENT_FILE $TMPDIR/fragment.intervals
cp $TARGET_FILE $TMPDIR/target.intervals


#GATK resources
INDELS_1000G=#indels1000G
INDELS_GOLDSTD=#indelsGoldStd
DBSNP=#dbSnp

echo "`${NOW}`INFO $SCRIPT_CODE copying GATK resources to tmp directory..."
INDELS_1000G_FILENAME=`basename $INDELS_1000G`
INDELS_GOLDSTD_FILENAME=`basename $INDELS_GOLDSTD`
DBSNP_FILENAME=`basename $DBSNP`

echo "`${NOW}`INFO $SCRIPT_CODE $INDELS_1000G"
cp $INDELS_1000G $TMPDIR/$INDELS_1000G_FILENAME
cp $INDELS_1000G.idx $TMPDIR/$INDELS_1000G_FILENAME.idx

echo "`${NOW}`INFO $SCRIPT_CODE $INDELS_GOLDSTD"
cp $INDELS_GOLDSTD $TMPDIR/$INDELS_GOLDSTD_FILENAME
cp $INDELS_GOLDSTD.idx $TMPDIR/$INDELS_GOLDSTD_FILENAME.idx

echo "`${NOW}`INFO $SCRIPT_CODE $DBSNP"
cp $DBSNP $TMPDIR/$DBSNP_FILENAME
cp $DBSNP.idx $TMPDIR/$DBSNP_FILENAME.idx

echo "`${NOW}`INFO $SCRIPT_CODE copying chunk BAM and index file to tmp directory..."
cp $INPUT_BAM $TMPDIR/chunk.bam
cp $INPUT_BAM.bai $TMPDIR/chunk.bam.bai

echo "`${NOW}`INFO $SCRIPT_CODE copying reference fasta and indexto tmp directory..."
cp $REFERENCE_FASTA $TMPDIR/reference.fa
cp $REFERENCE_FASTA.fai $TMPDIR/reference.fa.fai
cp $REFRENCE_SEQ_DICT $TMPDIR/reference.dict



# make tmp folder for temporary java files
mkdir $TMPDIR/tmp

# create target intervals for IndelRealigner !!!THIS LINE NOT TO BE CHANGED!!!!USED FOR MUTECT!!!
# although the input BAM file contains only a subset of reads
# RealignerTargetCreator will still traverse all sequences
# in the input header. Therefor, we still have to supply 
# the FRAGMENT_FILE telling it to process only the chunk
# for which the input BAM contains reads.
echo "`${NOW}`INFO $SCRIPT_CODE creating GATK realignment targets"
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
  -T RealignerTargetCreator -I $TMPDIR/chunk.bam \
  -nt $RTC_DATA_THREADS \
  -R $TMPDIR/reference.fa \
  -known $INDELS_1000G_FILENAME \
  -known $INDELS_GOLDSTD_FILENAME \
  -o $TMPDIR/$SAMPLE.$FRAGMENT.RTC.intervals \
  -L $TMPDIR/target.intervals 

echo "`${NOW}`INFO $SCRIPT_CODE copying realignment targets to $ANALYSIS_DIR/realignment/"
cp $TMPDIR/$SAMPLE.$FRAGMENT.RTC.intervals $ANALYSIS_DIR/realignment/

#logging
STATUS=OK
if [[ ! -s $ANALYSIS_DIR/realignment/$SAMPLE.$FRAGMENT.RTC.intervals ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\trealignment_targets\t$STATUS" >> $RUN_LOG

# run IndelRealigner
echo "`${NOW}`INFO $SCRIPT_CODE running GATK  realignment..."

INTERVAL_ARG="-L $TMPDIR/fragment.intervals"

if [ "$INCLUDES_UNMAPPED" == "T" ]
then
	INTERVAL_ARG="$INTERVAL_ARG -L unmapped"
fi 

java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	-R $TMPDIR/reference.fa \
	-I $TMPDIR/chunk.bam \
	-targetIntervals $SAMPLE.$FRAGMENT.RTC.intervals \
	-known $INDELS_1000G_FILENAME \
	-known $INDELS_GOLDSTD_FILENAME \
	-o $SAMPLE.$FRAGMENT.realigned.bam \
	$INTERVAL_ARG

# index realigned file
echo "`${NOW}`INFO $SCRIPT_CODE indexing realigned BAM..."
samtools index ${SAMPLE}.${FRAGMENT}.realigned.bam

# run BaseRecalibrator
echo "`${NOW}`INFO $SCRIPT_CODE generating recalibration report for realigned BAM..."
java -Xmx$JAVA_XMX -XX:+UseSerialGC -Djava.io.tmpdir=$TMPDIR/tmp -jar $GATK_HOME/GenomeAnalysisTK.jar \
   -nct $RTC_DATA_THREADS \
   -T BaseRecalibrator \
   -I ${SAMPLE}.${FRAGMENT}.realigned.bam \
   -R $TMPDIR/reference.fa \
   -knownSites $DBSNP_FILENAME \
   -knownSites $INDELS_1000G_FILENAME \
   -knownSites $INDELS_GOLDSTD_FILENAME \
   -o ${SAMPLE}.${FRAGMENT}.realigned.recal_data.grp \
   -L $TMPDIR/target.intervals \
   -rf BadCigar

echo "`${NOW}`INFO $SCRIPT_CODE copying recalibration report to $ANALYSIS_DIR/recalibration/reports/pre/..."
cp $SAMPLE.$FRAGMENT.realigned.recal_data.grp $ANALYSIS_DIR/recalibration/reports/pre/


#logging
STATUS=OK
if [[ ! -s $ANALYSIS_DIR/recalibration/reports/pre/$SAMPLE.$FRAGMENT.realigned.recal_data.grp ]]
then
	STATUS=FAILED
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\tprerecalibration_report\t$STATUS" >> $RUN_LOG


#step 5: soft clipping primer/probe sequences
if [[ $PRIMER_COORD_BED != ""  ]]
then
	echo "`${NOW}`INFO $SCRIPT_CODE soft clipping primer/probe sequences from realigned reads..."
	java -jar -Xmx$JAVA_XMX $NXTGENUTILS_HOME/NxtGenUtils.jar ClipPrimerSequences \
		-i ${SAMPLE}.${FRAGMENT}.realigned.bam \
		-o ${SAMPLE}.${FRAGMENT}.realigned.clipped.bam \
		-p $PRIMER_COORD_BED \
		-s $PRIMER_COORD_OFFSET
	
	#copy unclipped bam to output
	cp ${SAMPLE}.${FRAGMENT}.realigned.bam $ANALYSIS_DIR/realignment/${SAMPLE}.${FRAGMENT}.realigned.unclipped.bam
	cp ${SAMPLE}.${FRAGMENT}.realigned.bam.bai $ANALYSIS_DIR/realignment/${SAMPLE}.${FRAGMENT}.realigned.unclipped.bam.bai

	#logging
	STATUS=OK
	if [[ ! -s $TMPDIR/$SAMPLE.$FRAGMENT.realigned.clipped.bam ]]
	then
		STATUS=FAILED
	fi

	echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\tclipped_bam\t$STATUS" >> $RUN_LOG

				
	#replacing unclipped with clipped BAM
	mv ${SAMPLE}.${FRAGMENT}.realigned.clipped.bam ${SAMPLE}.${FRAGMENT}.realigned.bam
	
	#indexing clipped BAM
	samtools index ${SAMPLE}.${FRAGMENT}.realigned.bam
fi

# check that input and realigned bam have the same number of reads
# and if so, copy realigned bam to destination folder

READ_COUNT_INPUT=`samtools flagstat $TMPDIR/chunk.bam | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
READ_COUNT_OUTPUT=`samtools flagstat $TMPDIR/$SAMPLE.$FRAGMENT.realigned.bam | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`
echo "Number of reads in input chunk file: $READ_COUNT_INPUT"
echo "Number of reads in realigned file: $READ_COUNT_OUTPUT"

echo "`${NOW}`INFO $SCRIPT_CODE copying realigned BAM to output directory $ANALYSIS_DIR/realignment..."
cp $SAMPLE.$FRAGMENT.realigned.bam $ANALYSIS_DIR/realignment/
cp $SAMPLE.$FRAGMENT.realigned.bam.bai $ANALYSIS_DIR/realignment/


if [[ $READ_COUNT_INPUT -ne $READ_COUNT_OUTPUT ]]; then
	echo "WARNING: $SCRIPT_CODE input bam and realined bam have different numbers of reads"	
	echo "$SAMPLE $SCRIPT_CODE input chunk bam and realined bam have different numbers of reads" >> $WARNING_LOG
	echo "$SAMPLE $SCRIPT_CODE Number of reads in input chunk file: $READ_COUNT_INPUT" >> $WARNING_LOG
	echo "$SAMPLE $SCRIPT_CODE Number of reads in realigned chunk file: $READ_COUNT_OUTPUT" >> $WARNING_LOG
fi

#logging and deleting intermediate file

if [[ -s $ANALYSIS_DIR/realignment/$SAMPLE.$FRAGMENT.realigned.bam ]] && [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]] ; then
	STATUS=OK
	echo "deleting realigned bam file $INPUT_BAM"
	rm $INPUT_BAM $INPUT_BAM.bai
elif [[ -s $ANALYSIS_DIR/realignment/$SAMPLE.$FRAGMENT.realigned.bam ]] && [[ $READ_COUNT_INPUT -ne $READ_COUNT_OUTPUT ]]; then
	STATUS=WARNING
	echo "WARNING: $SCRIPT_CODE keeping input file for checks $INPUT_BAM"
else
	STATUS=FAILED
	echo "ERROR: $SCRIPT_CODE realignment failed, keeping input file for rerun $INPUT_BAM"
fi

echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\t$FRAGMENT\trealigned_bam\t$STATUS" >> $RUN_LOG

if  [[ "$STATUS" == "FAILED" ]]; then 
	echo "`${NOW}` realignement failed, exiting"
	exit 1
fi

#run summary script
perl $SUMMARY_SCRIPT_PATH

echo "`${NOW}`INFO $SCRIPT_CODE done"


ls -al


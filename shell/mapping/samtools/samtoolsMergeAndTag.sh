#!/bin/bash

#
# script to merge BWA alignments for
# a set reads
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=3:mem=20gb:tmpspace=tmpSpacemb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q queue

module load java/javaVersion
module load picard/picardVersion
module load samtools/samtoolsVersion
module load R/rVersion
module load gatk/gatkVersion
module load bedtools/bedtoolsVersion

PICARD_VERSION=picardVersion

#now
NOW="date +%Y-%m-%d%t%T%t"

OUTPUT_BAM_NAME=outBamName
PATH_OUTPUT_DIR=pathOutputDir
OUTPUT_BAM_PREFIX=`basename $OUTPUT_BAM_NAME .bam`

#java maximum memory
JAVA_XMX=19g

PICARD_MAX_RECORDS_IN_RAM=2000000

TMP_PATH_OUT_BAM=$TMPDIR/$OUTPUT_BAM_NAME
OUT_BAM=$PATH_OUTPUT_DIR/$OUTPUT_BAM_NAME

READ_GROUP_INFO=readGroupInfo
HEADER=headerFile
FILE2RG_MAPPING=file2RgMapping
TMP_PATH_RG_HEADER=$TMPDIR/rg_head.txt

PATH_BAIT_AMPLICON_INTERVALS=baitIntervalsFile
PATH_TARGET_INTERVALS=targetIntervalsFile
PATH_NON_OVERLAPPING_INTERVALS=nonOvelappingIntervalsFile
PATH_RIBOSOMAL_RNA_INTERVALS=ribosomalRnaIntervalsFile
PATH_ANNOTATION_REFFLAT=annotationRefFlat
REFERENCE_SEQUENCE=referenceSequence

CALCULATE_METRIC=calculateMetric
MARK_DUPLICATES=markDuplicates
METRIC_LEVEL=metricLevel
RG_LEVEL=`echo $METRIC_LEVEL | grep RG`
L_LEVEL=`echo $METRIC_LEVEL | grep L`
S_LEVEL=`echo $METRIC_LEVEL | grep S`
MAKE_BW=makeBw


echo "`${NOW}`-------------------------------------------------------------------------------------------------------"
#merge header and read group info
echo "`${NOW}`copying read group info to temp space..."
cp $READ_GROUP_INFO $TMP_PATH_RG_HEADER

echo "`${NOW}`copying reference sequence file to temp space..."
cp $REFERENCE_SEQUENCE $TMPDIR/tmp_reference.fa
cp $REFERENCE_SEQUENCE.fai $TMPDIR/tmp_reference.fa.fai

REFERENCE_SEQUENCE_PATH=`dirname $REFERENCE_SEQUENCE`
REFERENCE_SEQUENCE_BASENAME=`basename $REFERENCE_SEQUENCE .fa`
cp $REFERENCE_SEQUENCE_PATH/$REFERENCE_SEQUENCE_BASENAME.dict $TMPDIR/tmp_reference.dict

# Copy BAM files to merge to tmp space.
# The destination file name will be the
# read group name as samtools merge uses
# the filename as the read group name
# when tagging reads during merging. 
echo "`${NOW}`copying input BAM files to $TMPDIR..."

# Iterate over records in file.
# The first column contains the path to the
# BAM file to merge. The second column contains
# the read group name to be used to tag
# reads in the merged BAM file. samtools
# will use the filename as the read group tag. 
# Therefore the destination filename has to be
# the same as the read group name.
TOTAL_NUM_READS=0
while read FILE
do
  
  #get path to input BAM
  SOURCE=`echo $FILE | perl -e '$line=<>; \
      @tokens=split(/\s/,$line); \
      print @tokens[0];'`

  #get read group name
  DEST=`echo $FILE | perl -e '$line=<>; \
      @tokens=split(/\s/,$line); \
      print @tokens[1];'`

  #temp BAM file name
  DEST=$TMPDIR/$DEST.bam

  echo "`${NOW}`$SOURCE"
  echo "`${NOW}`to $DEST"

  cp $SOURCE $DEST
  rm $SOURCE


  NUM_READS=`samtools view -c $DEST`
  TOTAL_NUM_READS=$(( $TOTAL_NUM_READS + $NUM_READS ))

done < $FILE2RG_MAPPING

echo "`${NOW}`total number of reads in all input BAM files: $TOTAL_NUM_READS"

# We first merge BAM files on library level as
# marking of duplicates is performed on library level.
# Given correct read group tagging Picard tools
# should take into account library information when
# marking duplicates. However, it seems this is not
# the case and information on duplicate marking
# with Picard in context of libraries information
# could not be found in the documentation or 
# elsewhere online.


#for each library in the read group in file
LIBRARY_MERGED_BAM=""
for LIBRARY in `cut -f4 $TMP_PATH_RG_HEADER | uniq`
do

	IN_BAM=""

	#remove LB tag
	LIBRARY=`echo $LIBRARY | perl -pe 's/LB://g'`
	echo "`${NOW}`Merging library $LIBRARY"

	# we have to use the .dupmark.bam file extension as
	# it will be the duplicate marked BAMs that will get
	# merged
	LIBRARY_MERGED_BAM="$LIBRARY_MERGED_BAM INPUT=$TMPDIR/$LIBRARY.dupmark.bam"

	#generate file containing the RG header entries for the read
	#groups belonging to this library
	echo -n "" > library_rg_header.sam

	#get read group/BAM file names
	for READ_GROUP in `grep -P "LB:$LIBRARY\t" $TMP_PATH_RG_HEADER | cut -f2`
	do

		#remove ID tag
		READ_GROUP=`echo $READ_GROUP | perl -pe 's/ID://g'`

		# add BAM files to input parameter
		# for samtools merge
		IN_BAM="$IN_BAM $TMPDIR/$READ_GROUP.bam"

		# add RG header entry to RG header SAM for library
		grep -P "ID:$READ_GROUP\t" $TMP_PATH_RG_HEADER >> library_rg_header.sam

		#echo "`${NOW}`DEBUG: header of input BAM file $READ_GROUP.bam:"
   		#samtools view -H $TMPDIR/$READ_GROUP.bam
		
		if [ "$RG_LEVEL" ]
		then

		       	echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

			echo "`${NOW}`collecting base call quality and alignment summary metrics on read group level..."
		        #alignment summary metrics, insert size metrics, quality score distribution, mean quality by cycle
			OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$TMPDIR/tmp_reference.fa ASSUME_SORTED=true VERBOSITY=WARNING"
			java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CollectMultipleMetrics.jar INPUT=$TMPDIR/$READ_GROUP.bam OUTPUT=$TMPDIR/$READ_GROUP.bam $OPTIONS PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle

			echo "`${NOW}`copying quality metrics to output directory..."
			for METRIC in alignment_summary_metrics quality_by_cycle_metrics quality_by_cycle.pdf quality_distribution_metrics quality_distribution.pdf
			do
			        echo "`${NOW}`${OUTPUT_BAM_PREFIX}_${READ_GROUP}.$METRIC"
				cp $TMPDIR/$READ_GROUP.bam.$METRIC $PATH_OUTPUT_DIR/${OUTPUT_BAM_PREFIX}_${READ_GROUP}.$METRIC
				chmod 660 $PATH_OUTPUT_DIR/${OUTPUT_BAM_PREFIX}_${READ_GROUP}.$METRIC
			done;

			echo "`${NOW}`-------------------------------------------------------------------------------------------------------"
		fi 

	done;

	#generate file containing the header dictionary of the last
	#BAM file that was seen above with the RG header entries
	#appended
	samtools view -H $TMPDIR/$READ_GROUP.bam | grep -vP '^@RG' > library_header_dict.sam
	
	if [ "$MARK_DUPLICATES" = "TRUE" ]
	then
		#add Program header entry for Picard MarkDuplicates 
		echo -e "@PG\tID:PicardMarkDuplicates\tPN:MarkDuplicates\tVN:$PICARD_VERSION" >> library_header_dict.sam	
	fi

	#add Program header entry for Picard CleanSam
	echo -e "@PG\tID:PicardCleanSam\tPN:CleanSam\tVN:$PICARD_VERSION" >> library_header_dict.sam	
	
	#concatenate header dictionary and RG records
	cat library_rg_header.sam library_header_dict.sam > library_header.sam	
	  

	#get BAM count 
	IN_BAM_COUNT=`echo $IN_BAM | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`
	echo "`${NOW}`merging bam files $IN_BAM"

	#merge BAM files
	echo "`${NOW}`merging $IN_BAM_COUNT BAM files for library $LIBRARY..."

	# if there is more then one BAM file
	# belonging to the library
	if [ $IN_BAM_COUNT -ge 2 ]
    	then		
		    
    		echo "`${NOW}`samtools merge -rh $TMPDIR/library_header.sam $LIBRARY.bam $IN_BAM"    
#    		samtools merge -rh $TMPDIR/library_header.sam $LIBRARY.bam $IN_BAM
    		samtools merge -h $TMPDIR/library_header.sam $LIBRARY.bam $IN_BAM
       
	fi

	# if there is only one BAM file
	# belonging to the library
	if [ $IN_BAM_COUNT -eq 1 ]
    	then

		echo "`${NOW}`only $IN_BAM_COUNT BAM file for library $LIBRARY. Nothing to merge."

	        #we still need to run samtools merge 
                #to add the read group tag:

		#create an empty BAM file
		samtools view -bH $IN_BAM > dummy.bam

		#run merge to add read group tag
#		samtools merge -rh $TMPDIR/library_header.sam $LIBRARY.bam $IN_BAM dummy.bam
		samtools merge -h $TMPDIR/library_header.sam $LIBRARY.bam $IN_BAM dummy.bam

	fi

	if [ $IN_BAM_COUNT -eq 0 ]
    	then
		
		echo "`${NOW}`WARNING! $IN_BAM_COUNT BAM files for library $LIBRARY."

	fi

	#echo "`${NOW}`DEBUG: header of merged BAM file:"
	#samtools view -H $LIBRARY.bam

	for BAM in $IN_BAM
	do
	    rm $BAM
	done

	echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

	if [ "$MARK_DUPLICATES" = "TRUE" ]
	then
	
		echo "`${NOW}`marking duplicates on library level..."
		OPTIONS="METRICS_FILE=$TMPDIR/$LIBRARY.dupmark.stats ASSUME_SORTED=true TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT PROGRAM_RECORD_ID=null VERBOSITY=WARNING"
		java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/MarkDuplicates.jar INPUT=$TMPDIR/$LIBRARY.bam OUTPUT=$TMPDIR/$LIBRARY.dupmark.bam $OPTIONS 

		echo "`${NOW}`removing unmarked BAM..."
		rm $TMPDIR/$LIBRARY.bam

		echo "`${NOW}`copying duplicate marking stats to $LIBRARY.dupmark.stats..."
		cp $TMPDIR/$LIBRARY.dupmark.stats $PATH_OUTPUT_DIR/${OUTPUT_BAM_PREFIX}_${LIBRARY}.dupmark.stats
		chmod 660 $PATH_OUTPUT_DIR/${OUTPUT_BAM_PREFIX}_${LIBRARY}.dupmark.stats

		echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

	else
		
		echo "`${NOW}`skipping duplicate marking..."

		#just rename BAM file
		mv $TMPDIR/$LIBRARY.bam $TMPDIR/$LIBRARY.dupmark.bam

	fi


	if [ "$L_LEVEL" ]
	then

	        #GC bias
		echo "`${NOW}`collecting GC bias metrics on library level..."
		OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=tmp_reference.fa ASSUME_SORTED=true VERBOSITY=WARNING"
		java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar	$PICARD_HOME/CollectGcBiasMetrics.jar INPUT=$TMPDIR/$LIBRARY.dupmark.bam OUTPUT=$TMPDIR/$LIBRARY.dupmark.bam.gcBias CHART_OUTPUT=$TMPDIR/$LIBRARY.dupmark.bam.gcBias.pdf SUMMARY_OUTPUT=$TMPDIR/$LIBRARY.dupmark.bam.gcBiasSummary $OPTIONS

		#insert size metrics
		echo "`${NOW}`collecting insert size metrics on library level..."
		java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar	$PICARD_HOME/CollectInsertSizeMetrics.jar INPUT=$TMPDIR/$LIBRARY.dupmark.bam OUTPUT=$TMPDIR/$LIBRARY.dupmark.bam.insert_size_metrics HISTOGRAM_FILE=$TMPDIR/$LIBRARY.dupmark.bam.insert_size_histogram.pdf METRIC_ACCUMULATION_LEVEL=ALL_READS $OPTIONS

		echo "`${NOW}`copying output metrics to output directory..."
		for METRIC in gcBias gcBias.pdf gcBiasSummary insert_size_metrics insert_size_histogram.pdf 
		do
		    echo "`${NOW}`${OUTPUT_BAM_PREFIX}_${LIBRARY}.$METRIC"
		    cp $TMPDIR/$LIBRARY.dupmark.bam.$METRIC $PATH_OUTPUT_DIR/${OUTPUT_BAM_PREFIX}_${LIBRARY}.$METRIC
		    chmod 660 $PATH_OUTPUT_DIR/${OUTPUT_BAM_PREFIX}_${LIBRARY}.$METRIC
		done;
		echo "`${NOW}`-------------------------------------------------------------------------------------------------------"
        fi        
		
done;

# count number of library BAMs to merge
IN_BAM_COUNT=`echo $LIBRARY_MERGED_BAM | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`

#do final merge...
echo "`${NOW}`merging $IN_BAM_COUNT library BAM files and replacing header..."

if [ $IN_BAM_COUNT -ge 2 ]
then
    	
	echo "`${NOW}`merging library BAM files and replacing header..."
	#echo "	java -jar -Xmx$JAVA_XMX $PICARD_HOME/MergeSamFiles.jar $LIBRARY_MERGED_BAM OUTPUT=$TMP_PATH_OUT_BAM.sorted.dupmark SORT_ORDER=coordinate USE_THREADING=true  VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR"
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/MergeSamFiles.jar $LIBRARY_MERGED_BAM OUTPUT=$TMP_PATH_OUT_BAM.sorted.dupmark SORT_ORDER=coordinate USE_THREADING=true  VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR VERBOSITY=WARNING
        
fi


if [ $IN_BAM_COUNT -eq 1 ]
then

	echo "`${NOW}`only $IN_BAM_COUNT library BAM file. Nothing to merge."

	#nothing to merge, just rename file
	LIBRARY_MERGED_BAM=`echo $LIBRARY_MERGED_BAM | perl -pe 's/INPUT=//'`
	mv $LIBRARY_MERGED_BAM $TMP_PATH_OUT_BAM.sorted.dupmark

fi

#echo "`${NOW}`DEBUG: header of sample BAM file $TMP_PATH_OUT_BAM.sorted.dupmark:"
#samtools view -H $TMP_PATH_OUT_BAM.sorted.dupmark

echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

echo "`${NOW}`cleaning BAM file..."
OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT VERBOSITY=WARNING"
java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CleanSam.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark OUTPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean $OPTIONS 

echo "`${NOW}`removing unclean BAM..."
rm $TMP_PATH_OUT_BAM.sorted.dupmark

echo "`${NOW}`indexing sorted merged BAM"
samtools index $TMP_PATH_OUT_BAM.sorted.dupmark.clean

TOTAL_NUM_READS_AFTER_MERGE=`samtools view -c $TMP_PATH_OUT_BAM.sorted.dupmark.clean`
echo "`${NOW}`total number of reads in sorted merged BAM $TOTAL_NUM_READS_AFTER_MERGE"
if [ $TOTAL_NUM_READS_AFTER_MERGE != $TOTAL_NUM_READS_AFTER_MERGE ]
then
    echo "`${NOW}`Number of reads before and after merging is not the same!!!!!!!!!!!!!!!"
    exit 1
fi

echo "`${NOW}`creating flagstat for sorted merged BAM"
samtools flagstat $TMP_PATH_OUT_BAM.sorted.dupmark.clean > $TMP_PATH_OUT_BAM.flagstat

echo "`${NOW}`copying merged, dupmarked and clean BAM to $OUT_BAM"
cp $TMP_PATH_OUT_BAM.sorted.dupmark.clean $OUT_BAM
chmod 660 $OUT_BAM

echo "`${NOW}`copying BAM index to $OUT_BAM.bai"
cp $TMP_PATH_OUT_BAM.sorted.dupmark.clean.bai $OUT_BAM.bai
chmod 660 $OUT_BAM.bai

echo "`${NOW}`copying flagstats to $OUT_BAM.flagstats..."
cp $TMP_PATH_OUT_BAM.flagstat $OUT_BAM.flagstat
chmod 660 $OUT_BAM.flagstat

echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

if [ "$S_LEVEL" ]
then

        #generate metrics
	echo "`${NOW}`collecting quality metrics from merged file"

        #alignment summary metrics, insert size metrics, quality score distribution, mean quality by cycle
	OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$TMPDIR/tmp_reference.fa ASSUME_SORTED=true VERBOSITY=WARNING"
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CollectMultipleMetrics.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean OUTPUT=$TMP_PATH_OUT_BAM $OPTIONS PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle

        #GC bias
	OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=tmp_reference.fa ASSUME_SORTED=true VERBOSITY=WARNING"
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CollectGcBiasMetrics.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean OUTPUT=$TMP_PATH_OUT_BAM.gcBias CHART_OUTPUT=$TMP_PATH_OUT_BAM.gcBias.pdf SUMMARY_OUTPUT=$TMP_PATH_OUT_BAM.gcBiasSummary $OPTIONS 

	echo "`${NOW}`copying quality metrics to output directory..."
	for METRIC in alignment_summary_metrics gcBias gcBias.pdf gcBiasSummary insert_size_metrics insert_size_histogram.pdf quality_by_cycle_metrics quality_by_cycle.pdf quality_distribution_metrics quality_distribution.pdf
	do
	    echo "`${NOW}`$OUT_BAM.$METRIC"
	    cp $TMP_PATH_OUT_BAM.$METRIC $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.$METRIC
	    chmod 660 $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.$METRIC
	done;

	echo "`${NOW}`-------------------------------------------------------------------------------------------------------"
fi

#generate application specific 
#metrics if applicable
if [ "$CALCULATE_METRIC" == "HS" ]
then

	echo "`${NOW}`collecting Hybrid Sequencing metrics from merged file"
	echo "`${NOW}`using bait coordinates  : $PATH_BAIT_AMPLICON_INTERVALS"
	echo "`${NOW}`using target coordinates: $PATH_TARGET_INTERVALS"
	echo "`${NOW}`using reference sequence file: $REFERENCE_SEQUENCE"

        echo "`${NOW}`copying bait interval file to temp space..."
	cp $PATH_BAIT_AMPLICON_INTERVALS $TMPDIR/tmp_bait.int

        echo "`${NOW}`copying target interval file to temp space..."
	cp $PATH_TARGET_INTERVALS $TMPDIR/tmp_target.int
	
        echo "`${NOW}`collecting metrics..."	
	OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$TMPDIR/tmp_reference.fa VERBOSITY=WARNING"
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CalculateHsMetrics.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean OUTPUT=$TMP_PATH_OUT_BAM.hybridMetrics PER_TARGET_COVERAGE=$TMP_PATH_OUT_BAM.perTargetCoverage BAIT_INTERVALS=$TMPDIR/tmp_bait.int TARGET_INTERVALS=$TMPDIR/tmp_target.int $OPTIONS

	echo "`${NOW}`copying metrics to $OUT_BAM.hybridMetrics..."
	cp $TMP_PATH_OUT_BAM.hybridMetrics $OUT_BAM.hybridMetrics
	chmod 660 $OUT_BAM.hybridMetrics

	echo "`${NOW}`copying per target coverage to $OUT_BAM.perTargetCoverage..."
	cp $TMP_PATH_OUT_BAM.perTargetCoverage $OUT_BAM.perTargetCoverage
	chmod 660 $OUT_BAM.perTargetCoverage

	echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

elif [ "$CALCULATE_METRIC" == "TP" ]
then

	echo "`${NOW}`collecting Targeted PCR metrics from merged file"
	echo "`${NOW}`using amplicon coordinates: $BAIT_AMPLICON_INTERVALS"
	echo "`${NOW}`using target coordinates  : $TARGET_INTERVALS"
	echo "`${NOW}`using reference sequence file: $REFERENCE_SEQUENCE"

	echo "`${NOW}`copying bait interval file to temp space..."
	cp $PATH_BAIT_AMPLICON_INTERVALS $TMPDIR/tmp_amplicon.intervals

	echo "`${NOW}`copying target interval file to temp space..."
	cp $PATH_TARGET_INTERVALS $TMPDIR/tmp_target.intervals


	#INFO: failes without meaningful error message if fasta file index (generated with samtools faidx)
	#is not present.
	echo "`${NOW}`collecting Picard Tartgeted PCR metrics..."	
	OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=tmp_reference.fa VERBOSITY=WARNING"
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CollectTargetedPcrMetrics.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean OUTPUT=$TMP_PATH_OUT_BAM.targetedPcrMetrics PER_TARGET_COVERAGE=$TMP_PATH_OUT_BAM.perTargetCoverage AMPLICON_INTERVALS=$TMPDIR/tmp_amplicon.intervals TARGET_INTERVALS=$TMPDIR/tmp_target.intervals $OPTIONS

	echo "`${NOW}`copying metrics to $OUT_BAM.targetedPcrMetrics..."
	cp $TMP_PATH_OUT_BAM.targetedPcrMetrics $OUT_BAM.targetedPcrMetrics
	chmod 660 $OUT_BAM.targetedPcrMetrics

	echo "`${NOW}`copying all per target coverage to $OUT_BAM.perTargetCoverage..."
	cp $TMP_PATH_OUT_BAM.perTargetCoverage $OUT_BAM.perTargetCoverage
	chmod 660 $OUT_BAM.perTargetCoverage

	# if we provided non-overlapping regions of the amplicon file
	# generate coverage statistics for it

	if [[ "$PATH_NON_OVERLAPPING_INTERVALS" != "none" ]]; then
		echo "`${NOW}`copying non-overlapping amplicon interval file to temp space..."
		cp $PATH_NON_OVERLAPPING_INTERVALS $TMPDIR/tmp_amplicon.non_overlapping.intervals

		echo "`${NOW}`collecting Picard Tartgeted PCR metrics for non-overlapping amplicon intervals..."	
		OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=tmp_reference.fa VERBOSITY=WARNING"
		java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CollectTargetedPcrMetrics.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean OUTPUT=$TMP_PATH_OUT_BAM.non_overlapping.targetedPcrMetrics PER_TARGET_COVERAGE=$TMP_PATH_OUT_BAM.non_overlapping.perTargetCoverage AMPLICON_INTERVALS=$TMPDIR/tmp_amplicon.non_overlapping.intervals TARGET_INTERVALS=$TMPDIR/tmp_amplicon.non_overlapping.intervals $OPTIONS

		echo "`${NOW}`copying non-overlapping per target coverage to $OUT_BAM.perTargetCoverage..."
		cp $TMP_PATH_OUT_BAM.non_overlapping.perTargetCoverage $OUT_BAM.non_overlapping.perTargetCoverage
		chmod 660 $OUT_BAM.non_overlapping.perTargetCoverage
	fi

	#change filename of input file as GATK only
	#accepts input BAM files with the .bam extension :-S
	mv $TMP_PATH_OUT_BAM.sorted.dupmark.clean $TMP_PATH_OUT_BAM.sorted.dupmark.clean.bam
	
	echo "`${NOW}`collecting GATK amplicon coverage metrics..."	
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $GATK_HOME/GenomeAnalysisTK.jar \
		-R tmp_reference.fa \
   		-T DepthOfCoverage \
   		--countType COUNT_FRAGMENTS \
   		--omitDepthOutputAtEachBase \
   		-o $TMP_PATH_OUT_BAM.depthOfCoverage \
   		-I $TMP_PATH_OUT_BAM.sorted.dupmark.clean.bam \
		-rf BadCigar \
   		-L $TMPDIR/tmp_amplicon.intervals

	for EXT in sample_cumulative_coverage_counts sample_cumulative_coverage_proportions sample_interval_statistics sample_interval_summary sample_statistics sample_summary
	do
		cp $TMP_PATH_OUT_BAM.depthOfCoverage.$EXT $OUT_BAM.perAmpliconCoverage.$EXT
		chmod 660 $OUT_BAM.perAmpliconCoverage.$EXT
	done;

	# if using non-overlapping intervals file, calculate gatk coverage for non-overlapping regions
 
	if [[ "$PATH_NON_OVERLAPPING_INTERVALS" != "none" ]]; then

		echo "`${NOW}`collecting GATK amplicon coverage metrics..."	
		java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $GATK_HOME/GenomeAnalysisTK.jar \
			-R tmp_reference.fa \
   			-T DepthOfCoverage \
   			--countType COUNT_FRAGMENTS \
   			--omitDepthOutputAtEachBase \
   			-o $TMP_PATH_OUT_BAM.non_overlapping.depthOfCoverage \
   			-I $TMP_PATH_OUT_BAM.sorted.dupmark.clean.bam \
			-rf BadCigar \
   			-L $TMPDIR/tmp_amplicon.non_overlapping.intervals

		for EXT in sample_cumulative_coverage_counts sample_cumulative_coverage_proportions sample_interval_statistics sample_interval_summary sample_statistics sample_summary
		do
			cp $TMP_PATH_OUT_BAM.non_overlapping.depthOfCoverage.$EXT $OUT_BAM.non_overlapping.perAmpliconCoverage.$EXT
			chmod 660 $OUT_BAM.non_overlapping.perAmpliconCoverage.$EXT
		done;

	fi


	echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

elif [ "$CALCULATE_METRIC" == "RS" ]
then

	echo "`${NOW}`collecting RNA-seq metrics from merged file"
        echo "`${NOW}`using ribosomal RNA coordinates: $PATH_RIBOSOMAL_RNA_INTERVALS"
        echo "`${NOW}`using annotation refFlat: $PATH_ANNOTATION_REFFLAT"

        echo "`${NOW}`copying ribosomal RNA coordinates file to temp space..."
	cp $PATH_RIBOSOMAL_RNA_INTERVALS $TMPDIR/ribosomal.int

        echo "`${NOW}`copying annotation refFlat file to temp space..."
	cp $PATH_ANNOTATION_REFFLAT $TMPDIR/annotation.refFlat

        echo "`${NOW}`collecting metrics..."	
	OPTIONS="TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE=$TMPDIR/tmp_reference.fa VERBOSITY=WARNING"
	java -XX:+UseSerialGC -Xmx$JAVA_XMX -jar $PICARD_HOME/CollectRnaSeqMetrics.jar INPUT=$TMP_PATH_OUT_BAM.sorted.dupmark.clean OUTPUT=$TMP_PATH_OUT_BAM.RnaSeqMetrics CHART_OUTPUT=$TMP_PATH_OUT_BAM.chartOutput RIBOSOMAL_INTERVALS=$TMPDIR/ribosomal.int REF_FLAT=$TMPDIR/annotation.refFlat STRAND_SPECIFICITY=NONE $OPTIONS

	echo "`${NOW}`copying metrics to $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.RnaSeqMetrics...."
	cp $TMP_PATH_OUT_BAM.RnaSeqMetrics $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.RnaSeqMetrics
	chmod 660 $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.RnaSeqMetrics

	echo "`${NOW}`copying plot of normalized position vs. coverage to $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.chartOutput..."
	cp $TMP_PATH_OUT_BAM.chartOutput $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.chartOutput
	chmod 660 $PATH_OUTPUT_DIR/$OUTPUT_BAM_PREFIX.chartOutput

	echo "`${NOW}`-------------------------------------------------------------------------------------------------------"

fi

if [ "$MAKE_BW" == "TRUE" ]
then

	cat $TMPDIR/tmp_reference.dict | perl -e 'while(<>) {if (/^\@SQ\tSN:(\d+|X|Y|MT)\tLN:(\d+)\t.*/) {$chrom = $1; $size = $2; $chrom = "M" if $chrom eq "MT"; print "chr$chrom\t$size\n"}}' > $TMPDIR/chrom.sizes

	echo "`${NOW}` generating bigWig file to see coverage profile in UCSC browser"

	echo "`${NOW}` removing duplicates and unmapped reads and renaming chromosome names to hg19" 

	samtools view -hF 1028 $TMP_PATH_OUT_BAM.sorted.dupmark.clean | perl -e 'while(<>) {if (/^\@SQ\tSN:(\d+|X|Y)/) { s/(.*\tSN:)(.*)/$1chr$2/; print } elsif (/^\@SQ\tSN:MT/) { s/(.*\tSN:)MT(.*)/$1chrM$2/; print } elsif (/^\@SQ/ && /SN:/) { next } elsif (/^@/) { print } else { @cols = split(/\t/); if ($cols[2] =~ /^\d+$|^X$|^Y$/) {$cols[2] = "chr$cols[2]"} elsif ($cols[2] =~ /^MT$/) {$cols[2] = "chrM"} else {next} $line = ""; foreach $col (@cols) { $line .= $col."\t"; } chop($line); print "$line" }}' | samtools view -bS - > $TMPDIR/$OUTPUT_BAM_PREFIX.filtered.renamed.bam
	samtools index $TMPDIR/$OUTPUT_BAM_PREFIX.filtered.renamed.bam
	samtools flagstat $TMPDIR/$OUTPUT_BAM_PREFIX.filtered.renamed.bam > $TMPDIR/$OUTPUT_BAM_PREFIX.filtered.renamed.flagstat

	TOTAL=`cat $TMPDIR/$OUTPUT_BAM_PREFIX.filtered.renamed.flagstat|grep total|cut -f 1 -d ' '`
	SCALE=`echo "scale=10;1000000/${TOTAL}" | bc`

	echo "`${NOW}` calculate coverage, multiply by scale factor $SCALE ($TOTAL reads)"
	genomeCoverageBed -bg -split -scale $SCALE -ibam $TMPDIR/$OUTPUT_BAM_PREFIX.filtered.renamed.bam -g $TMPDIR/chrom.sizes > $TMPDIR/$OUTPUT_BAM_PREFIX.bedGraph

	echo "`${NOW}` translate bedGraph to BigWig"
	/groupvol/cgi/software/ucsc/bedGraphToBigWig $TMPDIR/$OUTPUT_BAM_PREFIX.bedGraph $TMPDIR/chrom.sizes $TMPDIR/$OUTPUT_BAM_PREFIX.bw
	cp $TMPDIR/$OUTPUT_BAM_PREFIX.bw $PATH_OUTPUT_DIR

fi

#list files for debugging
ls -alh

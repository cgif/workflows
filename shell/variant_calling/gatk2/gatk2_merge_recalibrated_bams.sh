#!/bin/bash

## script to merge realigned and recalibrated BAM files

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=5gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi


# load required modules
module load picard/#picardVersion
module load samtools/#samtoolsVersion
module load java/#javaVersion

NOW="date +%Y-%m-%d%t%T%t"
JAVA_XMX=4G

SCRIPT_CODE="GATKMEBA"

# define variables
IN_BAM_REALIGNED="#realignedBam"
IN_BAM_RECALIBRATED="#recalibratedBam"
IN_BAM_REDUCED="#reducedBam"
INPUT_DIR_REALIGNED=#inputDirRealigned
INPUT_DIR_RECALIBRATED=#inputDirRecalibrated
OUTPUT_PREFIX=#outputPrefix
PATH_OUTPUT_DIR_REALIGNED=#pathOutputDirRealigned
PATH_OUTPUT_DIR_RECALIBRATED=#pathOutputDirRecalibrated
PATH_STRIP_INDEL_QUALS_SCRIPT=#pathStripIndelQualsScript
SAMPLE=#sample
RUN_LOG=#runLog

#do not merge $IN_BAM_REALIGNED

#for each set of BAM files
#we will only merge realigned and recalibrated 
#BAM files to save storage space
for IN_BAM in "$IN_BAM_RECALIBRATED"
do

	OUT_BAM_NAME=$OUTPUT_PREFIX
	PATH_OUTPUT_DIR=$PATH_OUTPUT_DIR_REALIGNED
	
	INPUT_DIR=$INPUT_DIR_REALIGNED
	if [[ $IN_BAM == *recalibrated* ]]
	then
		
		INPUT_DIR=$INPUT_DIR_RECALIBRATED
		PATH_OUTPUT_DIR=$PATH_OUTPUT_DIR_RECALIBRATED
		
		OUT_BAM_NAME=$OUT_BAM_NAME.recalibrated
		
		if [[ $IN_BAM == *reduced* ]]
		then
			OUT_BAM_NAME=$OUT_BAM_NAME.reduced
		fi
		
	fi
	
	OUT_BAM_NAME=$OUT_BAM_NAME.bam
	
	IN_BAM_COUNT=`echo $IN_BAM | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`
	TMP_PATH_OUT_BAM=$TMPDIR/tmp.bam


	#as we only merge the realigned-recalibrated 
	#BAM files we will save the merged file simply
	#under the sample name
	#OUT_BAM=$PATH_OUTPUT_DIR/$OUT_BAM_NAME
	OUT_BAM=$PATH_OUTPUT_DIR/$SAMPLE.bam

	READ_COUNT_INPUT=0
	READ_COUNT_OUTPUT=0

	# if there is more than
	# one input BAM file
	if [ $IN_BAM_COUNT -ge 2 ]
	then

    	# copy BAM files to be merged to temp space and get total number of reads in all bam files before merging
	    echo "`${NOW}`INFO $SCRIPT_CODE copying BAM files to $TMPDIR..."
	    TMP_IN_BAM=""
	    INPUT_READ_COUNT=0

	    for BAM in $IN_BAM
	    do		
	
			BAM_BASENAME=`basename $BAM`
			echo "`${NOW}`INFO $SCRIPT_CODE $BAM_BASENAME"
			cp $INPUT_DIR/$BAM_BASENAME $TMPDIR
			TMP_IN_BAM="$TMP_IN_BAM INPUT=$TMPDIR/$BAM_BASENAME"

			#get number of reads in the input BAM file
	        READ_COUNT=`samtools flagstat $TMPDIR/$BAM_BASENAME | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	        READ_COUNT_INPUT=$(($READ_COUNT_INPUT + $READ_COUNT))

	    done
	    
	    # merge BAM files with picard tools
	    # picard allows to merge unsorted BAM files and
	    # output a coordinate sorted BAM while samtools will output
	    # an unsorted BAM file if the input files are unsorted.
	    # (Increasing the number of reads in memory before spilling to disc
	    #  has no impact on the runtime! This was tested with 10M reads which
	    # requires ~8GB of RAM and 4 processors. The default is 500k reads.)	
			
	    echo "`${NOW}`INFO $SCRIPT_CODE merging and coordiante sorting BAM files..."
	    java -jar -Xmx$JAVA_XMX $PICARD_HOME/MergeSamFiles.jar $TMP_IN_BAM OUTPUT=$TMP_PATH_OUT_BAM SORT_ORDER=coordinate USE_THREADING=true  VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR
    
    	#remove indel quality scores from recalibrated BAM files
		if [[ $IN_BAM == *recalibrated* ]] ||
		   [[ $IN_BAM == *reduced* ]]
		then
	    	echo "`${NOW}`INFO $SCRIPT_CODE stripping BAM file of indel quality scores..."
    		$PATH_STRIP_INDEL_QUALS_SCRIPT -i $TMP_PATH_OUT_BAM -o $TMP_PATH_OUT_BAM.stripped
    		mv $TMP_PATH_OUT_BAM.stripped $TMP_PATH_OUT_BAM
        fi
    
    	# get number of reads in the output bam file
	    READ_COUNT_OUTPUT=`samtools flagstat $TMP_PATH_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
    
	    # index output BAM
		echo "`${NOW}`INFO $SCRIPT_CODE indexing merged BAM..."
		samtools index $TMP_PATH_OUT_BAM
	    
	    # copy merged BAM, index, flagstat report and index
		# to destination folder
		echo "`${NOW}`INFO $SCRIPT_CODE copying merged BAM to $OUT_BAM..."
		cp $TMP_PATH_OUT_BAM $OUT_BAM
		chmod 660 $OUT_BAM
	
		echo "`${NOW}`INFO $SCRIPT_CODE copying BAM index to $OUT_BAM.bai"
		cp $TMP_PATH_OUT_BAM.bai $OUT_BAM.bai
		chmod 660 $OUT_BAM.bai
    
    	#logging
    	$BAM=""
    	if [[ $IN_BAM == *recalibrated* ]]
    	then
    		$BAM=recalibrated_bam_merged
    	elif [[ $IN_BAM == *reduced* ]]
    	then
    		$BAM=reduced_bam_merged
    	fi
    	
		STATUS=FAILED
		if [[ ! -e $OUT_BAM ]]
		then
			STATUS=FAILED
		fi
		echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t$BAM\t$STATUS" >> $RUN_LOG

    
	fi

	# if there is only one input BAM
	# nothing to merge. Just copy BAM
	# file to temp space for sorting, indexing
	# flagstat and md5
	if [ $IN_BAM_COUNT -eq 1 ]
	then

	    echo "`${NOW}`INFO $SCRIPT_CODE only one input BAM file. Nothing to merge."
	    BAM_BASENAME=`basename $IN_BAM`
		IN_BAM=$INPUT_DIR/$BAM_BASENAME
   
    	#remove indel quality scores from recalibrated BAM files
		if [[ $IN_BAM == *recalibrated* ]] ||
		   [[ $IN_BAM == *reduced* ]]
		then
		
			echo "`${NOW}`INFO $SCRIPT_CODE copying input BAM to $TMPDIR..."
	    	cp $IN_BAM $TMPDIR

			# get number of reads in the input bam file
		    READ_COUNT_INPUT=`samtools flagstat $BAM_BASENAME | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	    			
	    	#strip indel qualities
	    	echo "`${NOW}`INFO $SCRIPT_CODE stripping BAM file of indel quality scores..."
    		$PATH_STRIP_INDEL_QUALS_SCRIPT -i $BAM_BASENAME -o $TMP_PATH_OUT_BAM.stripped
    		mv $TMP_PATH_OUT_BAM.stripped $TMP_PATH_OUT_BAM
  
  			# get number of reads in the output bam file
		    READ_COUNT_OUTPUT=`samtools flagstat $TMP_PATH_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	    	
  			# index output BAM
			echo "`${NOW}`INFO $SCRIPT_CODE indexing stripped BAM..."
			samtools index $TMP_PATH_OUT_BAM
         
         	echo "`${NOW}`INFO $SCRIPT_CODE copying merged BAM to $OUT_BAM..."
			cp $TMP_PATH_OUT_BAM $OUT_BAM
			chmod 660 $OUT_BAM
         
         	echo "`${NOW}`INFO $SCRIPT_CODE copying BAM index to $OUT_BAM.bai"
			cp $TMP_PATH_OUT_BAM.bai $OUT_BAM.bai
			chmod 660 $OUT_BAM.bai
         
        else
        
        	# get number of reads in the input bam file
		    READ_COUNT_INPUT=`samtools flagstat $IN_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	    	
        	echo "`${NOW}`INFO $SCRIPT_CODE copying input BAM to $OUT_BAM..."
	    	cp $IN_BAM $OUT_BAM
	    	chmod 660 $OUT_BAM

		    READ_COUNT_OUTPUT=`samtools flagstat $OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
		
		 	echo "`${NOW}`INFO $SCRIPT_CODE copying input BAM index to $OUT_BAM.bai..."	
	    	cp $IN_BAM.bai $OUT_BAM.bai
	    	chmod 660 $OUT_BAM.bai
    
        fi
    	   
    	#logging
    	BAM=""
    	if [[ $IN_BAM == *recalibrated* ]]
    	then
    		BAM=recalibrated_bam_merged
    	fi
    	
    	if [[ $IN_BAM == *reduced* ]]
    	then
    		BAM=reduced_bam_merged
    	fi
    	
		STATUS=OK
		if [[ ! -e $OUT_BAM ]]
		then
			STATUS=FAILED
		fi
		echo -e "`${NOW}`$SCRIPT_CODE\t$SAMPLE\tall\t$BAM\t$STATUS" >> $RUN_LOG
    	   
	fi

	echo "`${NOW}`INFO $SCRIPT_CODE input  BAM(s) read count: $READ_COUNT_INPUT"
	echo "`${NOW}`INFO $SCRIPT_CODE output BAM    read count: $READ_COUNT_OUTPUT"
	
	#...if yes, delete input BAM files
	if [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]] &&
	   [[ -e $OUT_BAM ]]
	then
	
		echo "`${NOW}`INFO $SCRIPT_CODE deleting intermediate BAM files..."
		rm $IN_BAM

	#...if no, keep input BAM files for re-run
	else
	
		echo "`${NOW}`WARN $SCRIPT_CODE Output BAM does not contain the same number of reads as the input BAM file(s)!"
		echo "`${NOW}`WARN $SCRIPT_CODE Keeping intermediate BAM files for re-run..."  		 
	
	fi

done;


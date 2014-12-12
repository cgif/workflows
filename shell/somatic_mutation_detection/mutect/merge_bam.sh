#!/bin/bash

## script to merge realigned and recalibrated BAM files

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=5gb:tmpspace=#tmpSpaceGbgb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load required modules
module load samtools/#samtoolsVersion

NOW="date +%Y-%m-%d%t%T%t"

# define variables
IN_BAM="#recalibratedBam"
INPUT_PATH=#inputDirRecalibrated
OUTPUT_PATH=#pathOutputDirRecalibrated
SAMPLE=#sample

IN_BAM_COUNT=`echo $IN_BAM | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`
TMP_IN_BAM=""
TMP_OUT_BAM=$TMPDIR/tmp.bam
OUT_BAM=$OUTPUT_PATH/$SAMPLE.bam

READ_COUNT_INPUT=0
READ_COUNT_OUTPUT=0

# if there is more than one input BAM file
if [ $IN_BAM_COUNT -ge 2 ]; then

	echo "`${NOW}`copying BAM files to $TMPDIR..."
	for BAM in $IN_BAM; do 

		# copy BAM files to be merged to temp space
		echo "`${NOW}`$BAM"
		cp $INPUT_PATH/$BAM $TMPDIR

		#get number of reads in the input BAM file
		READ_COUNT=`samtools flagstat $TMPDIR/$BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
		READ_COUNT_INPUT=$(($READ_COUNT_INPUT + $READ_COUNT))

    		#remove indel quality scores from recalibrated BAM files
		echo "`${NOW}`stripping BAM file of indel quality scores..."
		samtools view -h $TMPDIR/$BAM | perl -pe 's/B[ID]:Z:.*?\t//g' | samtools view -bS - >  $TMPDIR/$BAM.stripped
		mv $TMPDIR/$BAM.stripped $TMPDIR/$BAM

		TMP_IN_BAM="$TMP_IN_BAM $TMPDIR/$BAM"
	
	done
	ls -lh
	du -hs $TMPDIR

	#merge bam files
	echo "`${NOW}`samtools merge $TMP_OUT_BAM $TMP_IN_BAM"
        samtools merge $TMP_OUT_BAM $TMP_IN_BAM
	ls -lh
	du -hs $TMPDIR
   
	# get number of reads in the output bam file
	READ_COUNT_OUTPUT=`samtools flagstat $TMP_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	
	echo "`${NOW}`input  BAM(s) read count: $READ_COUNT_INPUT"
	echo "`${NOW}`output BAM    read count: $READ_COUNT_OUTPUT"
	    
	# copy merged BAM and index to destination folder
	if [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]]; then

		echo "`${NOW}`copying merged BAM to $OUT_BAM..."
		cp $TMP_OUT_BAM $OUT_BAM
		chmod 660 $OUT_BAM
         
		echo "`${NOW}`copying BAM index to $OUT_BAM.bai"
		samtools index $TMP_OUT_BAM
		cp $TMP_OUT_BAM.bai $OUT_BAM.bai
		chmod 660 $OUT_BAM.bai

	else

		echo "`${NOW}`Output BAM does not contain the same number of reads as the input BAM file(s)!"
		exit 1
 	
	fi

fi

#if there is only one input BAM nothing to merge. 
if [ $IN_BAM_COUNT -eq 1 ]; then
   
	#remove indel quality scores from recalibrated BAM files
	echo "`${NOW}`copying input BAM to $TMPDIR..."
	cp $INPUT_PATH/$IN_BAM $TMP_OUT_BAM

	# get number of reads in the input bam file
	READ_COUNT_INPUT=`samtools flagstat $TMP_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	

    	#remove indel quality scores from recalibrated BAM files
	echo "`${NOW}`stripping BAM file of indel quality scores..."
	samtools view -h $TMP_OUT_BAM | perl -pe 's/B[ID]:Z:.*?\t//g' | samtools view -bS - >  $TMP_OUT_BAM.stripped
	mv $TMP_OUT_BAM.stripped $TMP_OUT_BAM
   
	# get number of reads in the output bam file
	READ_COUNT_OUTPUT=`samtools flagstat $TMP_OUT_BAM | head -n 1 | perl -e 'while(<>){ if(/(\d*?)\s\+\s(\d*?)\s/) { $retval=$1+$2; print "$retval\n"; }  }'`	

	echo "`${NOW}`input  BAM(s) read count: $READ_COUNT_INPUT"
	echo "`${NOW}`output BAM    read count: $READ_COUNT_OUTPUT"
	    
	# copy merged BAM and index to destination folder

	if [[ $READ_COUNT_INPUT -eq $READ_COUNT_OUTPUT ]]; then

		echo "`${NOW}`copying merged BAM to $OUT_BAM..."
		cp $TMP_OUT_BAM $OUT_BAM
		chmod 660 $OUT_BAM
         
		echo "`${NOW}`copying BAM index to $OUT_BAM.bai"
		samtools index $TMP_OUT_BAM
		cp $TMP_OUT_BAM.bai $OUT_BAM.bai
		chmod 660 $OUT_BAM.bai

	else

		echo "`${NOW}`Output BAM does not contain the same number of reads as the input BAM file(s)!"
		exit 1
 	
	fi

fi	


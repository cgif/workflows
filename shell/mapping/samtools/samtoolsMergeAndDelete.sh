#!/bin/bash

#
# script to merge a set of BAM files
# the output BAM file is coordinate sorted
# BAM index, flagstat report and md5 sum
# are generated for the output BAM

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=2:mem=4800mb:tmpspace=100gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi


module load samtools/0.1.19
module load picard/1.85

#now
NOW="date +%Y-%m-%d%t%T%t"

GROUP_VOL_CGI=/groupvol/cgi

OUTPUT_PREFIX=outputPrefix
PATH_OUTPUT_DIR=pathOutputDir
INPUT_DIR=inputDir
IN_BAM=inBam
IN_BAM_COUNT=`echo $IN_BAM | perl -e '$in=<>; @tokens=split(/\s/,$in); $count=@tokens; print $count;'`

TMP_PATH_OUT_BAM=$TMPDIR/tmp.bam
OUT_BAM=$PATH_OUTPUT_DIR/$OUTPUT_PREFIX.sorted.bam

# if there is more than
# one input BAM file
if [ $IN_BAM_COUNT -ge 2 ]
then

    # copy BAM files to be merged to temp space and get total number of reads in all bam files before merging
    echo "`${NOW}`copying BAM files to $TMPDIR..."
    TMP_IN_BAM=""
    INPUT_READ_COUNT=0

    for BAM in $IN_BAM
    do
    		
		BAM_BASENAME=`basename $BAM`
		echo "`${NOW}`$BAM_BASENAME"
		cp $INPUT_DIR/$BAM $TMPDIR
		TMP_IN_BAM="$TMP_IN_BAM INPUT=$TMPDIR/$BAM_BASENAME"

        READ_COUNT=`samtools flagstat $TMPDIR/$BAM_BASENAME | head -n 1 | cut -f1 -d ' '`	
        INPUT_READ_COUNT=$(($INPUT_READ_COUNT + $READ_COUNT))

        if [ "$READ_COUNT" -eq 0 ] || [ -z $READ_COUNT ]
        then
            echo "`$NOW`file $BAM_BASENAME contains no reads." 
            exit 1
        fi

    done
    
    echo "`$NOW`total number of reads in the entire set of intermediate bam files before merging: $INPUT_READ_COUNT"
    
    # merge BAM files with picard tools
    # picard allows to merge unsorted BAM files and
    # output a coordinate sorted BAM while samtools will output
    # an unsorted BAM file if the input files are unsorted.
    # (Increasing the number of reads in memory before spilling to disc
    #  has no impact on the runtime! This was tested with 10M reads which
    # requires ~8GB of RAM and 4 processors. The default is 500k reads.)	
		
    echo "`${NOW}`merging and coordiante sorting BAM files..."
    java -jar -Xmx4000m $PICARD_HOME/MergeSamFiles.jar $TMP_IN_BAM OUTPUT=$TMP_PATH_OUT_BAM SORT_ORDER=coordinate USE_THREADING=true  VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR

fi

# if there is only one input BAM
# nothing to merge. Just copy BAM
# file to temp space for sorting, indexing
# flagstat and md5
if [ $IN_BAM_COUNT -eq 1 ]
then

#   get number of reads in the input bam file
    IN_BAM=`echo $IN_BAM | perl -pe 's/\s//g'`

    INPUT_READ_COUNT=`samtools flagstat $INPUT_DIR/$IN_BAM | head -n 1 | cut -f1 -d ' '`	
    echo "`$NOW`number of reads in the input bam file: $INPUT_READ_COUNT" 

    if [ "$INPUT_READ_COUNT" -eq 0 ] || [ -z $INPUT_READ_COUNT ]
    then
        echo "`$NOW`file $INPUT_DIR/$IN_BAM contains no reads." 
        exit 1
    fi

    echo "`${NOW}`only one input BAM file. Nothing to merge."
    cp $INPUT_DIR/$IN_BAM $TMP_PATH_OUT_BAM

    echo "`${NOW}`coordinate sorting input BAM..." 	
    #sort BAM file
    samtools sort -m 5000000 $TMP_PATH_OUT_BAM $TMP_PATH_OUT_BAM.sorted
    
    #replace unsorted BAM with sorted BAM
    mv $TMP_PATH_OUT_BAM.sorted.bam $TMP_PATH_OUT_BAM

fi

# clean BAM (soft-clip an alignment that hangs off the end 
# of its reference sequence and set MAPQ to 0 if read is unmapped)
echo "`${NOW}`cleaning merged BAM..."
java -jar -Xmx4000m $PICARD_HOME/CleanSam.jar INPUT=$TMP_PATH_OUT_BAM OUTPUT=$TMP_PATH_OUT_BAM.clean VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR
mv $TMP_PATH_OUT_BAM.clean $TMP_PATH_OUT_BAM

# index output BAM
echo "`${NOW}`indexing merged BAM..."
samtools index $TMP_PATH_OUT_BAM

# create flagstat report
echo "`${NOW}`creating flagstat..."
samtools flagstat $TMP_PATH_OUT_BAM > $TMP_PATH_OUT_BAM.flagstat

echo "`${NOW}`creating md5 checksum..."
md5sum $TMP_PATH_OUT_BAM > $TMP_PATH_OUT_BAM.md5

# copy merged BAM, index, flagstat report and index
# to destination folder
echo "`${NOW}`copying merged BAM to $OUT_BAM..."
cp $TMP_PATH_OUT_BAM $OUT_BAM
chmod 640 $OUT_BAM

echo "`${NOW}`copying BAM index to $OUT_BAM.bai"
cp $TMP_PATH_OUT_BAM.bai $OUT_BAM.bai
chmod 640 $OUT_BAM.bai

echo "`${NOW}`copying flagstat to $OUT_BAM.flagstat"
cp $TMP_PATH_OUT_BAM.flagstat $OUT_BAM.flagstat
chmod 640 $OUT_BAM.flagstat

echo "`${NOW}`copying md5sum to $OUT_BAM.md5"
cp $TMP_PATH_OUT_BAM.md5 $OUT_BAM.md5
chmod 640 $OUT_BAM.md5

#make sure that output has correct number of reads...
OUTPUT_READ_COUNT=`cat $TMP_PATH_OUT_BAM.flagstat | head -n 1 | cut -f 1 -d ' '`

echo "`${NOW}`input reads (read1 + read2): $INPUT_READ_COUNT"
echo "`${NOW}`reads in output: $OUTPUT_READ_COUNT"

#...if yes, delete input BAM files
if [ $OUTPUT_READ_COUNT -eq $INPUT_READ_COUNT ]
then
	echo "`${NOW}`deleting intermediate BAM files..."
	rm $IN_BAM

#...if no, keep input BAM files for re-run
else
	echo "`${NOW}`WARNING!!! Output BAM does not contain the same number of reads as the input files!"
        echo "`${NOW}`keeping intermediate BAM files for re-run..."  		 
fi


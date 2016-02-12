#!/bin/bash

#script to run rins for identification of non-human sequences 

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=#threadsPerRun
#PBS -l mem=50g
#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#load modules
module load blast+/#blastVersion
module load blat/#blatVersion
module load trinity/#trinityVersion
module load bowtie/#bowtieVersion

BASEDIR=#baseDir
CONFIG=#configPath
RESULT=#resultPath

BOWTIE_INDEX_HOST=#bowtieIndexHost
BLASTN_INDEX_HOST=#blastIndexHost
BLASTN_INDEX_VIRAL=#blastIndexViral
FASTA_VIRAL=#fastaViral
PATH_READS_DIRECTORY=#readsDir
PATTERN_READ1=#patternRead1
PATTERN_READ2=#patternRead2

#copy reference sequences to tmp space
BOWTIE_INDEX_HOST_BASENAME=`basename $BOWTIE_INDEX_HOST`
cp $BOWTIE_INDEX_HOST.1.ebwt $TMPDIR/
cp $BOWTIE_INDEX_HOST.2.ebwt $TMPDIR/
cp $BOWTIE_INDEX_HOST.3.ebwt $TMPDIR/
cp $BOWTIE_INDEX_HOST.4.ebwt $TMPDIR/
cp $BOWTIE_INDEX_HOST.rev.1.ebwt $TMPDIR/
cp $BOWTIE_INDEX_HOST.rev.2.ebwt $TMPDIR/

BLASTN_INDEX_HOST_BASENAME=`basename $BLASTN_INDEX_HOST`
cp $BLASTN_INDEX_HOST.nhr $TMPDIR/
cp $BLASTN_INDEX_HOST.nin $TMPDIR/
cp $BLASTN_INDEX_HOST.nsq $TMPDIR/

BLASTN_INDEX_VIRAL_BASENAME=`basename $BLASTN_INDEX_VIRAL`
cp $BLASTN_INDEX_VIRAL.nhr $TMPDIR/
cp $BLASTN_INDEX_VIRAL.nin $TMPDIR/
cp $BLASTN_INDEX_VIRAL.nsq $TMPDIR/

FASTA_VIRAL_BASENAME=`basename $FASTA_VIRAL`
cp $FASTA_VIRAL $TMPDIR/


TMP_PATH_READS_FASTQ_READ1=""
TMP_PATH_READS_FASTQ_READ2=""
#for each read1 fastq file 
for FASTQ_READ1 in `ls --color=never $PATH_READS_DIRECTORY/*.f*q* | grep $PATTERN_READ1`
do
 
        FASTQ_READ1=`basename $FASTQ_READ1`

        #find read2 mate file
        FASTQ_READ2=""
	for FASTQ in `ls --color=never $PATH_READS_DIRECTORY/*.f*q* | grep $PATTERN_READ2`
	do	

	        FASTQ=`basename $FASTQ`
    		FASTQ_REPLACE=`echo $FASTQ | perl -pe "s/$PATTERN_READ2/$PATTERN_READ1/"`

    		if [ "$FASTQ_REPLACE" = "$FASTQ_READ1" ]; 
    		then
		        FASTQ_READ2=$FASTQ     
	    	fi

	done;
             
        #check if mate file found and the number of lines in mate files is the same
	if [ -z $FASTQ_READ2 ]
	then
	        echo "No mate file found for $FASTQ_READ1. Skipped."   		
	else

                COUNT_LINES_READ1=`gzip -d -c $PATH_READS_DIRECTORY/$FASTQ_READ1 | wc -l | cut -f 1 -d ' '`
                COUNT_LINES_READ2=`gzip -d -c $PATH_READS_DIRECTORY/$FASTQ_READ2 | wc -l | cut -f 1 -d ' '`

                if [ $COUNT_LINES_READ1 -eq $COUNT_LINES_READ2 ]
                then

                        #replace empty and sequence lines (together with corresponding quality lines) created by cutadapt
                        #unzip and copy fastqs to tmp space
                        echo "`${NOW}`copying reads to temporary scratch space..."

                        #generate string that will replace empty and short sequence lines in fastq files
                        READ_LENGTH=`gzip -d -c $PATH_READS_DIRECTORY/$FASTQ_READ1 | head -n 100 | awk '{if(NR%4==2) print length($1)}' | sort -n| uniq | tail -n 1`
                        STRING=$(for i in `eval echo {1..$READ_LENGTH}`;do printf "%s" "N";done;)

                        FASTQ_READ1_NO_EXT=`basename $FASTQ_READ1 .gz`
	                echo "`${NOW}`$FASTQ_READ1_NO_EXT"
                        gzip -c -d $PATH_READS_DIRECTORY/$FASTQ_READ1 | sed "s/^$/$STRING/g" > $TMPDIR/$FASTQ_READ1_NO_EXT
	                TMP_PATH_READS_FASTQ_READ1="$TMP_PATH_READS_FASTQ_READ1 $TMPDIR/$FASTQ_READ1_NO_EXT"

                        FASTQ_READ2_NO_EXT=`basename $FASTQ_READ2 .gz`
	                echo "`${NOW}`$FASTQ_READ2_NO_EXT"
                        gzip -c -d $PATH_READS_DIRECTORY/$FASTQ_READ2 | sed "s/^$/$STRING/g" > $TMPDIR/$FASTQ_READ2_NO_EXT
	                TMP_PATH_READS_FASTQ_READ2="$TMP_PATH_READS_FASTQ_READ2 $TMPDIR/$FASTQ_READ2_NO_EXT"

                else
                        echo "Unequal number of lines in the mate files. Skipped." 
                fi

	fi

done

#remove trailing white space
TMP_PATH_READS_FASTQ_READ1=`echo $TMP_PATH_READS_FASTQ_READ1 | perl -pe 's/^ //'`
TMP_PATH_READS_FASTQ_READ2=`echo $TMP_PATH_READS_FASTQ_READ2 | perl -pe 's/^ //'`

sed -i -e "s/#tmpBowtieIndexHost/${BOWTIE_INDEX_HOST_BASENAME//\//\\/}/" $CONFIG
sed -i -e "s/#tmpBlastIndexHost/${BLASTN_INDEX_HOST_BASENAME//\//\\/}/" $CONFIG
sed -i -e "s/#tmpBlastIndexViral/${BLASTN_INDEX_VIRAL_BASENAME//\//\\/}/" $CONFIG
sed -i -e "s/#tmpFastaViral/${FASTA_VIRAL_BASENAME//\//\\/}/" $CONFIG
sed -i -e "s/#tmpRead1/${TMP_PATH_READS_FASTQ_READ1//\//\\/}/" $CONFIG
sed -i -e "s/#tmpRead2/${TMP_PATH_READS_FASTQ_READ2//\//\\/}/" $CONFIG


perl $BASEDIR/rins.pl -c $CONFIG -o $RESULT

ls -l $TMPDIR 
RESULT_DIR=`dirname $RESULT`
cp $TMPDIR/*leftlane* $RESULT_DIR
cp $TMPDIR/*rightlane* $RESULT_DIR
cp $TMPDIR/*human_contig* $RESULT_DIR
cp $TMPDIR/*non_human* $RESULT_DIR
cp $TMPDIR/*clean_blastn* $RESULT_DIR
cp $TMPDIR/*Trinity* $RESULT_DIR

chmod 660 $RESULT_DIR/*

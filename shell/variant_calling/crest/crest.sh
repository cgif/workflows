#!/bin/bash

## script to run CREST

#PBS -l walltime=walltimeHours:00:00
#PBS -l select=1:ncpus=1:mem=2gb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

# load modules
module load samtools/0.1.18
module load crest/1.0
module load blat/34

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

# define variables
REFERENCE_FASTA=referenceFasta
REFERENCE_2BIT=reference2bit
INPUT_BAM=inputBam
RESULTS_FOLDER=resultsFolder
READ_LENGTH=readLength
MIN_SC_READS=minScReads
SENSITIVE=sensitive
PORT=serverPort
DEPLOYMENT_SERVER=deploymentServer
DEPLOYMENT_PATH=deploymentPath

BAM_NAME=`basename $INPUT_BAM .bam`

# step 1: get soft-clipping positions
echo "`${NOW}`getting soft-clipping positions"
echo "`${NOW}`files $RESULTS_FOLDER/$BAM_NAME.bam.cover and $RESULTS_FOLDER/$BAM_NAME.bam.sclip.txt will be generated"

extractSClip.pl -i $INPUT_BAM --ref_genome $REFERENCE_FASTA -o $RESULTS_FOLDER

# step 2: run SV detection script
echo "`${NOW}`running SV detection script"
echo "`${NOW}`file $RESULTS_FOLDER/$BAM_NAME.bam.predSV.txt will be generated"

gfServer start localhost $PORT $REFERENCE_2BIT &

STAT=0
while [ $STAT -lt 2 ]
do
    sleep 20
    STAT=`ps -elf|grep "gfServer start localhost $PORT"|grep " S "|wc -l`
done 

# Run SV detection

if [[ "$SENSITIVE" == "yes" ]]; then
 	CREST.pl -f $RESULTS_FOLDER/${BAM_NAME}.bam.cover -d $INPUT_BAM --ref_genome $REFERENCE_FASTA -t $REFERENCE_2BIT --blatserver localhost --blatport $PORT -l $READ_LENGTH -o $RESULTS_FOLDER --min_sclip_reads $MIN_SC_READS --sensitive 
else
	CREST.pl -f $RESULTS_FOLDER/${BAM_NAME}.bam.cover -d $INPUT_BAM --ref_genome $REFERENCE_FASTA -t $REFERENCE_2BIT --blatserver localhost --blatport $PORT -l $READ_LENGTH -o $RESULTS_FOLDER --min_sclip_reads $MIN_SC_READS
fi

# Kill gfserver
gfServer stop localhost $PORT

# step 3: create html file to display alignment at breakpoint
echo "`${NOW}`creating breakpoint html file"
echo "`${NOW}`file $RESULTS_FOLDER/$BAM_NAME.bam.predSV.html will be generated"

bam2html.pl -d $INPUT_BAM -i $RESULTS_FOLDER/${BAM_NAME}.bam.predSV.txt --ref_genome $REFERENCE_FASTA -o $RESULTS_FOLDER/${BAM_NAME}.predSV.html

# change permissions of the results files
chmod 660 $RESULTS_FOLDER/*

#deply html file to the server
scp $RESULTS_FOLDER/${BAM_NAME}.predSV.html  $DEPLOYMENT_SERVER:$DEPLOYMENT_PATH/${BAM_NAME}.predSV.html
ssh $DEPLOYMENT_SERVER chmod -R 664 $DEPLOYMENT_PATH/${BAM_NAME}.predSV.html

echo "`${NOW}` done"


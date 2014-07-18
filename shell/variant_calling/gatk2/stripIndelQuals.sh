#!/bin/bash

#strips insertion (BI) and deletion (BD) qualities from 
#GATK recalibrated BAM files

#now
NOW="date +%Y-%m-%d%t%T%t"


USAGE="USAGE: stripIndelQuals.sh -i <input_bam> -o <output_bam>"

#parse command line args
while getopts "i:o:" OPTION;
do

    case "$OPTION" in
	i) INPUT_PATH="$OPTARG";;
	o) OUTPUT_PATH="$OPTARG";;
	h) echo $USAGE;;
	[?]) echo $USAGE;;

esac
done

#check if required args are set
if [[ -z $INPUT_PATH ]]
then

        #...if not print usage and exit
        echo "`$NOW`ERROR: required argument input BAM (-i) missing."
        echo $USAGE
        exit 1
fi

if [[ -z $OUTPUT_PATH ]]
then

        #...if not print usage and exit
        echo "`$NOW`ERROR: required argument output BAM (-o) missing."
        echo $USAGE
        exit 1
fi


#check if input BAM exists
if [[ ! -e $INPUT_PATH ]]
then

        #...if not throw error and exit
        echo "`$NOW`ERROR: input file does not exist: $INPUT_PATH"
        exit 1
fi

#check if output BAM exists
if [[ -e $OUTPUT_PATH ]]
then

        #...if yes throw error and exit
        echo "`$NOW`ERROR: output file does exist already. Refusing to overwrite."
        exit 1
fi

#decompress BAM, replace BD/BI annotations in SAM and recompress 
samtools view -h $INPUT_PATH | perl -pe 's/B[ID]:Z:.*?\t//g' | samtools view -bS - > $OUTPUT_PATH


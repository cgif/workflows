#!/bin/bash

#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=1gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

INPUT_DIR=#input_dir
OUTPUT_FILE=#output_file

NOW="date +%Y-%m-%d%t%T%t"

printf "SAMPLE\tFREQ\tLK_DIFF\n" > $OUTPUT_FILE

for SAMPLE in `ls --color=never $INPUT_DIR|grep -v multisample`; do

	FREQ=`cut -f 7 $INPUT_DIR/${SAMPLE}/${SAMPLE}.selfSM|grep -v FREE`
	FREQ=`echo "scale=2; $FREQ/1" | bc`
	LK1=`cut -f 8 $INPUT_DIR/${SAMPLE}/${SAMPLE}.selfSM|grep -v FREE`
	LK0=`cut -f 9 $INPUT_DIR/${SAMPLE}/${SAMPLE}.selfSM|grep -v FREE`
	LK_DIFF=`echo "scale=2; $LK1 - $LK0" | bc`

	printf "$SAMPLE\t$FREQ\t$LK_DIFF\n" >> $OUTPUT_FILE

done



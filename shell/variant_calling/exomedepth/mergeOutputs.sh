#!/bin/bash

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

# load modules

NOW="date +%Y-%m-%d%t%T%t"

# define variables
PED_FILE=#pedFile
MS_RESULTS_DIR=#resultsFolder
MAKE_BED=#makeBed
PROJECT=#project
SUMMARY_SCRIPT_PATH=#summaryScriptPath
COMBINE_CALLS_PATH=#combineCallsPath
EXONS_PATH=#exonsPath

RESULTS_DIR=`dirname $MS_RESULTS_DIR`
MERGED_CNVS=$MS_RESULTS_DIR/${PROJECT}.cnvs.tsv
MERGED_SUMMARY=$MS_RESULTS_DIR/${PROJECT}.cnvs.summary.tsv
COMBINED_CALLS=$MS_RESULTS_DIR/${PROJECT}.cnvs.exons.tsv

echo "`$NOW`merging output for all samples"
#get headers
HEAD_SAMPLE=`grep -v '^#' $PED_FILE | cut -f 2 | head -n 1`
head -n 1 $RESULTS_DIR/$HEAD_SAMPLE/${HEAD_SAMPLE}.cnv.calls.all.tsv > $MERGED_CNVS
head -n 1 $RESULTS_DIR/$HEAD_SAMPLE/${HEAD_SAMPLE}.cnv.calls.all.summary.tsv > $MERGED_SUMMARY
chmod 660 $MERGED_CNVS
chmod 660 $MERGED_SUMMARY

for SAMPLE in `grep -v '^#' $PED_FILE | cut -f 2`; do

	tail -n +2 $RESULTS_DIR/$SAMPLE/${SAMPLE}.cnv.calls.all.tsv >> $MERGED_CNVS
	tail -n +2 $RESULTS_DIR/$SAMPLE/${SAMPLE}.cnv.calls.autosomes.summary.tsv >> $MERGED_SUMMARY
done

## now add X_chromosome summaries (visually easier to compare between samples)

for SAMPLE in `grep -v '^#' $PED_FILE | cut -f 2`; do

	tail -n +2 $RESULTS_DIR/$SAMPLE/${SAMPLE}.cnv.calls.X_chromosome.summary.tsv >> $MERGED_SUMMARY
done

echo "`$NOW`make multisample bed file"

MERGED_BED=$MS_RESULTS_DIR/${PROJECT}.cnvs.bed

$MAKE_BED $MERGED_CNVS $MERGED_BED

perl $COMBINE_CALLS_PATH $MERGED_CNVS $EXONS_PATH > ${COMBINED_CALLS}.unsorted

head -n 1 ${COMBINED_CALLS}.unsorted > $COMBINED_CALLS
tail -n +2 ${COMBINED_CALLS}.unsorted | sort -n -k4,4 -k5,5 >> $COMBINED_CALLS
sed -i '/^$/d' $COMBINED_CALLS

rm ${COMBINED_CALLS}.unsorted


chmod 660 $MS_RESULTS_DIR/*

#run summary script
perl $SUMMARY_SCRIPT_PATH

ls -al

echo "`${NOW}`done"


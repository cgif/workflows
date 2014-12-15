#!/bin/bash

## script to merge CREST outputs

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=2gb

#PBS -m ea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#now
NOW="date +%Y-%m-%d%t%T%t"

CREST_CHUNKS_TXT="#crestChunkTxt"
CREST_CHUNKS_HTML="#crestChunkHTML"
RESULTS_DIR=#resultsPath
SAMPLE=`basename $RESULTS_DIR`

DEPLOYMENT_SERVER=#deploymentServer
SUMMARY_DEPLOYMENT=#summaryDeployment

cat $CREST_CHUNKS_TXT > $RESULTS_DIR/$SAMPLE.predSV.txt
cat $CREST_CHUNKS_HTML > $RESULTS_DIR/$SAMPLE.predSV.html

scp -r $RESULTS_DIR/$SAMPLE.predSV.html $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT/ > /dev/null 2>&1
ssh $DEPLOYMENT_SERVER "chmod -R 664 $SUMMARY_DEPLOYMENT/$SAMPLE.predSV.html" > /dev/null 2>&1

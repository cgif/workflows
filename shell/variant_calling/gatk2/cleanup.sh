#!/bin/bash

## script to clean up directories after GATK run
## only run folders containing job scripts and logs are left in analysis folder for each sample
## 

NOW="date +%Y-%m-%d%t%T%t"

ANALYSIS_DIR=#analysisDir
RESULTS_DIR=#resultsDir

echo "analysis folder to be cleaned: $ANALYSIS_DIR"
echo "results folder to be cleaned: $RESULTS_DIR"

echo "`$NOW`cleaning analysis folder"

echo "`$NOW`removing all chunks folders"
for CHUNKS_FOLDER in `find $ANALYSIS_DIR -name chunks`
do
    echo "`$NOW`removing $CHUNKS_FOLDER"
    rm -R $CHUNKS_FOLDER
done;

echo "`$NOW`removing all realignement folders"
for REALIGNMENT_FOLDER in `find $ANALYSIS_DIR -name realignment`
do
    echo "`$NOW`removing $REALIGNMENT_FOLDER"
    rm -R $REALIGNMENT_FOLDER
done;

echo "`$NOW`removing all recalibration folders"
for RECALIBRATION_FOLDER in `find $ANALYSIS_DIR -name recalibration`
do
    echo "`$NOW`removing $RECALIBRATION_FOLDER"
    rm -R $RECALIBRATION_FOLDER
done;

echo "`$NOW`removing unifiedgenotyper folder from multisample"
rm -R $ANALYSIS_DIR/multisample/unifiedgenotyper

echo "`$NOW`cleaning results folder"
echo "`$NOW`removing all realignement folders"
for REALIGNMENT_FOLDER in `find $RESULTS_DIR -name realignment`
do
    echo "`$NOW`removing $REALIGNMENT_FOLDER"
    rm -R $REALIGNMENT_FOLDER
done;

echo "`$NOW`finished removing files"




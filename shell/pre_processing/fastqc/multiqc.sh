#!/bin/bash
#
# script to run multiqc.sh
#

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=500mb

#PBS -m ea
#PBS -M igf@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load multiqc
source activate multiqc

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

#set up script
FASTQ_DIR=#pathReadsFastq
RUNS_DIR=#pathRunsDir
REPORT_DIR=#pathReportsDir
MS_REPORT_DIR=#pathMSReportsDir
DEPLOYMENT_SERVER=#deploymentServer
SUMMARY_DEPLOYMENT=#summaryDeployment

/apps/anaconda/2.4.1/envs/multiqc/bin/multiqc $REPORT_DIR -n $MS_REPORT_DIR/multiqc.html

chmod -Rf 775 $MS_REPORT_DIR/*
scp -r $MS_REPORT_DIR/* $DEPLOYMENT_SERVER:$SUMMARY_DEPLOYMENT/
ssh $DEPLOYMENT_SERVER chmod -R 775 $SUMMARY_DEPLOYMENT/*


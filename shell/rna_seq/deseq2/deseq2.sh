#!/bin/bash

#
# runs DESeq2
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=1
#PBS -l mem=10gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe


module load R/#rVersion

RESULT_DIR=#resultsDir
R_SCRIPT=#rScript 
#export GFF_PATH=#gffPath
#export GO_PATH=#goPath
#DEPLOYMENT_SERVER=#deploymentServer
#SUMMARY_DEPLOYMENT=#summaryDeployment
ANALYSYS=`basename $R_SCRIPT`

#run R script
R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log


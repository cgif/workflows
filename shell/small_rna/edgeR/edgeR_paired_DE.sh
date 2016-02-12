#!/bin/bash

#
# runs edgeR on paired samples
#

#PBS -l walltime=#walltimeHours:00:00
#PBS -l ncpus=1
#PBS -l mem=2gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

module load R/3.1.0

RESULT_DIR=#resultsDir
R_SCRIPT=#rScript 

R CMD BATCH --no-save --no-restore $R_SCRIPT ${R_SCRIPT}.log

chmod 660 $RESULT_DIR/*

#!/bin/bash

## script to run mutect

#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=1gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

NOW="date +%Y-%m-%d%t%T%t"

# define variables
RESULTS_DIR=#results_dir
RAW_STATS_PATH="#rawStatsFiles"

SAMPLE=`basename $RESULTS_DIR`
RAW_STATS_MERGED=$RESULTS_DIR/$SAMPLE.stats

cat $RAW_STATS_PATH|grep -P '^contig\tposition\t'|head -n 1 > $RAW_STATS_MERGED

cat $RAW_STATS_PATH|grep -vP '^#'|grep -vP '^contig\tposition\t' >> $RAW_STATS_MERGED
 
grep -v REJECT $RAW_STATS_MERGED > $RAW_STATS_MERGED.keep






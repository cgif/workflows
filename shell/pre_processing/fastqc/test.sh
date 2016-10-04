#!/bin/bash

#
# generates, configures and submits scripts to run FastQC on
# a single fastq file, a set of fastq files in a directory or
# a CGI project directory
# test

USAGE="./qfastqc.usage"


#configuration
PATH_PLACE_HOLDER=forwardSlash

#now
NOW="date +%Y-%m-%d%t%T%t"

#today
TODAY=`date +%Y-%m-%d`

BASEDIR=`dirname $0`
GROUP_VOL_CGI=/groupvol/pgm/home/cgi

#parse command line args
IS_PROJECT_DIR=false
OUTPUT_PATH=$GROUP_VOL/results/fastqc   
DEPLOY=false

while getopts "i:o:hpd" option; do
    case "$option" in
	i) INPUT_PATH="$OPTARG";;
	p) IS_PROJET_DIR="true";;
	o) OUTPUT_PATH=$OPTARG";;
	d) DEPLOY="true";;
        h) cat $USAGE;;
	[?]) cat $USAGE;;
esac
done

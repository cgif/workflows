#!/bin/bash

## script to run FATHHM-MKL

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

module load python/2.7.3 
module load tabix/0.2.6

VCF=$1
OUT=$2
DB=/groupvol/cgi/bin/fathmm-MKL/fathmm-MKL_Current.tab.gz

cp $VCF $TMPDIR/tmp.vcf
cp $DB $TMPDIR/db.tab.gz
cp $DB.tbi $TMPDIR/db.tab.gz.tbi

/groupvol/cgi/bin/fathmm-MKL/fathmm-MKL.py $TMPDIR/tmp.vcf $TMPDIR/out.txt $TMPDIR/db.tab.gz

cp $TMPDIR/out.txt $OUT

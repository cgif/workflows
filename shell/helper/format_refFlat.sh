#!/bin/bash

#script for formating GFF annotation file into refFlat

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

cp $GFF $TMPDIR/tmp.gff

perl /groupvol/cgi/src/helper/format_refFlat.pl $TMPDIR

cp $TMPDIR/tmp.refFlat $RF
chmod 770 $RF
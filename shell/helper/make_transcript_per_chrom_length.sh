#!/bin/bash

#script for calculating transcript lengths per 
#chromosome from RefFlat file
#output file is used for mergeand tag plots

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M cgi@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

cp $RF $TMPDIR/tmp.rf

perl /home/asoskins/workflows/shell/helper/make_transcript_per_chrom_length.pl $TMPDIR

cp $TMPDIR/tmp.lf $LF
chmod 770 $LF

#!/bin/bash

#script for calculating gene lengths from GFF file

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M igf@imperial.ac.uk
#PBS -m bea
#PBS -j oe

#PBS -q pqcgi

cp $GFF $TMPDIR/tmp.gff

perl /project/tgu/src/helper/make_gene_length.pl $TMPDIR

cp $TMPDIR/tmp.length $LF
chmod 770 $LF

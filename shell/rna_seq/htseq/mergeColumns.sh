#!/bin/bash

#
# make counts table
#
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=1gb:tmpspace=1gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

htseq_files='#htseqFiles'
htseq_samples='#htseqSamples'
htseq_table=#htseqTable

names_file=`echo $htseq_files | perl -e '$string=<>; @cols=split(/\s/,$string); print @cols[-1];'`
cut -f1 $names_file > $TMPDIR/tmp.table

for column in $htseq_files
do

        if [ -f $column ]
	then

                cut -f 2 $column > $TMPDIR/tmp.column
		paste $TMPDIR/tmp.table $TMPDIR/tmp.column > $TMPDIR/tmp.table2
		cp $TMPDIR/tmp.table2 $TMPDIR/tmp.table

	fi

done

echo "$htseq_samples" > $TMPDIR/tmp.htseq.counts
cat $TMPDIR/tmp.table >> $TMPDIR/tmp.htseq.counts
cp $TMPDIR/tmp.htseq.counts $htseq_table
chmod 0660 $htseq_table
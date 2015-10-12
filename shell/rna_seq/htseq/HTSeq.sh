#!/bin/bash

#
# runs HTSeq
#
#PBS -l walltime=#walltimeHours:00:00
#PBS -l select=1:ncpus=1:mem=20gb:tmpspace=50gb

#PBS -m bea
#PBS -M cgi@imperial.ac.uk
#PBS -j oe

#PBS -q pqcgi

module load python/#pythonVersion
module load samtools/#samtoolsVersion

bam_file=#bamFile
htseq_counts=#htseqCounts
gff_file=#gffFile
ss_library=#strand
feature=#feature

cp $bam_file $TMPDIR/tmp.bam
cp $bam_file.bai $TMPDIR/tmp.bam.bai
cp $gff_file $TMPDIR/tmp.gff

echo "Sorting bam by name"
samtools sort -n $TMPDIR/tmp.bam nameSrt
echo "Fixing information for mates"
samtools fixmate $TMPDIR/nameSrt.bam $TMPDIR/nameSrt.fixmate.bam 
echo "Translating bam to sam" 
samtools view -h $TMPDIR/nameSrt.fixmate.bam > $TMPDIR/nameSrt.fixmate.sam
echo "Calculating counts with HTSeq python -m HTSeq.scripts.count -f sam -r name -s $ss_library -a 10 -t exon -i $feature -m union $TMPDIR/nameSrt.fixmate.sam $TMPDIR/tmp.gff"
python -m HTSeq.scripts.count -f sam -r name -s $ss_library -a 10 -t exon -i $feature -m union $TMPDIR/nameSrt.fixmate.sam $TMPDIR/tmp.gff > $TMPDIR/tmp.HTSeq.counts

grep -vP "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique" $TMPDIR/tmp.HTSeq.counts > $htseq_counts
grep -P "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique" $TMPDIR/tmp.HTSeq.counts > $htseq_counts.discard

htseq_dir=`dirname $htseq_counts`
cp $TMPDIR/nameSrt.fixmate.bam $htseq_dir/nameSrt.fixmate.bam

chmod 0660 $htseq_dir/*

ls -l 








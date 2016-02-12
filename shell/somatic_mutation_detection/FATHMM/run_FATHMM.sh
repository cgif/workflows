#!/bin/bash

## script to run FATHHM-MKL

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

module load python/2.7.3 
module load tabix/0.2.6

VCF=/groupvol/cgi/results/$PROJECT/oncotator/$DATE/$PROJECT
cp $VCF.annotated.vcf $TMPDIR/annotated.vcf

TODAY=`date +%Y-%m-%d`
RESULTS=/groupvol/cgi/results/$PROJECT/FATHMM/$TODAY
mkdir -p $RESULTS
chmod -R 770 $RESULTS

DB=/groupvol/cgi/bin/fathmm-MKL/fathmm-MKL_Current.tab.gz
cp $DB $TMPDIR/db.tab.gz
cp $DB.tbi $TMPDIR/db.tab.gz.tbi

#filter lost of function (LOF), non-synonymous and silent variants
#filter common variance (minor allele frequency in 1000 genomes or ESP above 1%)
cat $TMPDIR/annotated.vcf | grep -vE '^#' | perl -e 'while(<>) {@data=split(/\;/,$_); $noncoding=1; $uncommon=1; foreach $key (@data) {if ($key =~ /^variant_classification=(Frame_Shift_Del|Frame_Shift_Ins|De_novo_Start_OutOfFrame|Start_Codon_SNP|Splice_Site|Nonsense_Mutation|Nonstop_Mutation|In_Frame_Del|In_Frame_Ins|Missense_Mutation|Silent)/) {$noncoding--} if (/1000Genome_AF=(.*)/) {$uncommon-- if $1 > 0.01} if (/ESP_MAF=(.*)/) {@AF = split(/\|/, $1); foreach $AF (@AF) {@AF_all = split(/,/, $AF); $uncommon-- if $AF_all[2] > 1;}}} if ($noncoding == 1 && $uncommon == 1) {print}}' > $TMPDIR/noncoding.uncommon.annotated.vcf

#filter indels and put chromosome, coordinate, reference and mutation alleles into coma-separated file
cat $TMPDIR/noncoding.uncommon.annotated.vcf | perl -e 'while(<>) {@data=split(/\;/,$_); $f_snv=0; foreach $key (@data) {$f_snv++ if $key =~ /genome_change=.*del|genome_change=.*ins/} print unless $f_snv;}' | perl -e 'while(<>) {@data=split(/\t/,$_); $ref = $1 if $data[3] =~ /^(\w).*/; $mut = $1 if $data[4] =~ /^(\w).*/; print "$data[0],$data[1],$ref,$mut\n"}'  > $TMPDIR/noncoding.uncommon.snv.csv

/groupvol/cgi/bin/fathmm-MKL/fathmm-MKL.py $TMPDIR/noncoding.uncommon.snv.csv $TMPDIR/noncoding.uncommon.snv.FATHMM_scores.txt $TMPDIR/db.tab.gz

#uncommon noncoding SNVs with FATHMM score above 0.5 and all four annotation tracks available
cat $TMPDIR/noncoding.uncommon.snv.FATHMM_scores.txt | grep -vE '#|No Prediction Found' | awk '($5 > 0.5 && $6 == "ALL") {print}' | cut -f 1,2 > $TMPDIR/noncoding.uncommon.snv.FATHMM_scores.selected.txt

cp $TMPDIR/noncoding.uncommon.annotated.vcf $RESULTS/$PROJECT.noncoding.uncommon.annotated.vcf
cp $TMPDIR/noncoding.uncommon.snv.csv $RESULTS/$PROJECT.noncoding.uncommon.snv.csv
cp $TMPDIR/noncoding.uncommon.snv.FATHMM_scores.txt $RESULTS/$PROJECT.noncoding.uncommon.snv.FATHMM_scores.txt
cp $TMPDIR/noncoding.uncommon.snv.FATHMM_scores.selected.txt $RESULTS/$PROJECT.noncoding.uncommon.snv.FATHMM_scores.selected.txt


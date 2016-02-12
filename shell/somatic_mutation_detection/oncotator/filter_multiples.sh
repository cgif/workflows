#!/bin/bash

#script to find genes with multipe LOF & nonsyn mutations across all patients

#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=5gb

#PBS -M igf@imperial.ac.uk
#PBS -m ea
#PBS -j oe

#PBS -q pqcgi

#RESULTS_DIR=/groupvol/cgi/results/johnson_glioma/oncotator/2014-11-27
#PROJECT=johnson_glioma

RESULTS_DIR=#resultsDir
PROJECT=#project

VCF=$RESULTS_DIR/$PROJECT.annotated.vcf
cp $VCF $TMPDIR/vcf.tmp

echo -n "" > $TMPDIR/list.tmp

for GENE in `cat $TMPDIR/vcf.tmp|grep -v '^#'|grep -E ';variant_classification=(Frame_Shift_Del|Frame_Shift_Ins|De_novo_Start_OutOfFrame|Start_Codon_SNP|Splice_Site|Nonsense_Mutation|Nonstop_Mutation|In_Frame_Del|In_Frame_Ins|Missense_Mutation)'|perl -e 'while(<>) {$f_common=0; @data=split(/\;/,$_); foreach $key(@data) {if ($key =~ /^1000Genome_AF=(.*)/) {$f_common++ if $1 > 0.01} if ($key =~ /^ESP_MAF=(.*)/) {@AF = split(/\|/, $1); foreach $AF (@AF) {@AF_all = split(/,/, $AF); $f_common++ if $AF_all[2] > 1;}}} print unless $f_common}'|cut -f 8|perl -e 'while(<>) {@data=split(/\;/,$_); foreach $key (@data) {print "$1\n" if $key =~ /^gene=(.*)/}}'|sort|uniq`; do

	cat $TMPDIR/vcf.tmp|grep -v '^#'|grep -E "gene=$GENE\;"|grep -E ';variant_classification=(Frame_Shift_Del|Frame_Shift_Ins|De_novo_Start_OutOfFrame|Start_Codon_SNP|Splice_Site|Nonsense_Mutation|Nonstop_Mutation|In_Frame_Del|In_Frame_Ins|Missense_Mutation)'|perl -e 'while(<>) {$f_common=0; $variants=""; @data=split(/\;/,$_); foreach $key(@data) {if ($key =~ /^1000Genome_AF=(.*)/) {$f_common++ if $1 > 0.01} if ($key =~ /^ESP_MAF=(.*)/) {@AF = split(/\|/, $1); foreach $AF (@AF) {@AF_all = split(/,/, $AF); $f_common++ if $AF_all[2] > 1;}} if ($key =~ /^set=(.*)/) {$variants = $1}} $variants =~ s/-/\n/g; print "$variants\n" unless $f_common}' > $TMPDIR/gene_mut.tmp

	COUNT_MUT=`cat $TMPDIR/gene_mut.tmp|wc -l`

	COUNT_PAT=`cat $TMPDIR/gene_mut.tmp|sort|uniq|wc -l`
		
	if [ $COUNT_MUT -ge 2 ]; then

		GENE_ID=`cat $TMPDIR/vcf.tmp|grep -v '^#'|grep -E "gene=$GENE\;"|perl -e 'while(<>) {@data=split(/\;/,$_); foreach $key(@data) {print "$1" if $key =~ /.*=(ENSG\d+)/; last if $key =~ /.*=ENSG\d+/;} last;}'`

		BIOTYPE=`cat $TMPDIR/vcf.tmp|grep -v '^#'|grep -E "gene=$GENE\;"|perl -e 'while(<>) {@data=split(/\;/,$_); foreach $key(@data) {print "$1" if $key =~ /^gene_type=(.*)/} last;}'`
		
		printf "$GENE\t$GENE_ID\t$BIOTYPE\t$COUNT_PAT\t$COUNT_MUT\n" >> $TMPDIR/list.tmp

	fi

done

printf "Gene_name\tGene_ID\tGene_type\tPatients\tMutations\n" > $TMPDIR/list_sorted.tmp
sort -n -r -k 4 $TMPDIR/list.tmp >> $TMPDIR/list_sorted.tmp
cp $TMPDIR/list_sorted.tmp $RESULTS_DIR/$PROJECT.multiples.xls



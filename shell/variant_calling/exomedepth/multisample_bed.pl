#!/usr/bin/perl -w
use strict;

my $calls = $ARGV[0];
my $out = $ARGV[1];

open(CNV, $calls) or die "Cannot open cnv input file: $!"; 
open(OUT, ">$out") or die "Cannot open output file: $!";;

print OUT "browser hide all\n";
print OUT "browser pack decipher hgmd clinvar iscaViewDetail knownGene refGene ensGene dgvPlus\n";
print OUT "browser position chr15:25,333,120-25,333,929\n";
print OUT 'track name="ExomeDepth CNV" description="CNVs identified by ExomeDepth" visibility="pack" itemRgb="On" db="hg19"';
print OUT "\n";

my @data;
my $sample;
my $chrom;
my $start;
my $end;
my $colour;

while(<CNV>){
	chomp();
	next if /^sample/;
	@data = split(/\t/,$_);
	$sample = $data[0];
	$chrom = $data[5];
	next if $chrom =~ /GL/;
	$chrom = "chr"."$chrom";
	$start = $data[6];
	$end = $data[7];
	$colour="255,0,0" if $data[2] eq "deletion";
	$colour="0,0,255" if $data[2] eq "duplication";
	print OUT "$chrom\t$start\t$end\t$sample\t0\t+\t$start\t$end\t$colour\n";
}

close(CNV) or die "Cannot close input file: $!";
close(OUT) or die "Cannot close output file: $!";



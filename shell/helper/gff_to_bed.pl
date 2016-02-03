#!/usr/bin/perl -w

#script that converts Ensembl gff into bed file with gene coordinates, Ensembl Gene IDs, gene names and gene biotypes

use strict;
use warnings;

if ($#ARGV != 1 ) {
	print "usage: perl gff_to_bed.pl input_GFF_file output_bed_file\n";
	exit;
}

my $gff = $ARGV[0];
my $gene_bed = $ARGV[1];

open(GFF, "$gff");
open(BED, ">$gene_bed");

my @starts = (); my @ends = ();
my @sorted_starts = (); my @sorted_ends = ();
my @data = ();
my $gene = ""; my $chrom = "";
my $start = ""; my $end = "";
my $last_gene = ""; my $last_chrom = "";

while (<GFF>){
	chomp();

	@data = split(/\t/,$_);
	next unless $data[2] eq "exon";
	$chrom = $data[0];
	$start = $data[3];
	$end = $data[4];
	$gene = "$1;$2;$3" if $data[8] =~ /gene_id \"(.*?)\"\;.*gene_name \"(.*?)\"\;.*gene_biotype \"(.*?)\"\;/;

	if ("$gene" ne "$last_gene" && $last_gene ne ""){
		@sorted_starts = sort {$a <=> $b} @starts;
		@sorted_ends = sort {$b <=> $a} @ends;
		print BED "$last_chrom\t$sorted_starts[0]\t$sorted_ends[0]\t$last_gene\n";
		@starts = (); @ends = ();
		@sorted_starts = (); @sorted_ends = ();
	}

	push(@starts, $start);
	push(@ends, $end);

	$last_gene = "$gene";
	$last_chrom = "$chrom";
}
@sorted_starts = sort {$a <=> $b} @starts;
@sorted_ends = sort {$b <=> $a} @ends;
print BED "$last_chrom\t$sorted_starts[0]\t$sorted_ends[0]\t$last_gene\n";


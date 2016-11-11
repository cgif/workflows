#!/usr/bin/perl -w
use strict;

#script that converts Mutect and SomaticIndelDetector vcf files into tsv file for input into Oncotator web server

if ($#ARGV != 1 ) {
	print "usage: perl vcf_to_oncotator_tsv.pl input_vcf_file output_tsv_file\n";
	exit;
}

my $vcf = $ARGV[0];
my $tsv = $ARGV[1];

open(VCF, "$vcf") or die "Unable to open input file: $!";
open(TSV, ">$tsv") or die "Unable to open output file: $!";

my @vcfRecord;
my $chrom;
my $vcfPos;
my $start;
my $end;
my $ref;
my $alt;
my $vcfRef;
my $vcfAlt;

print TSV "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele\n";

while (my $line = <VCF> ){
	$line =~ /^#/ and next;
	chomp($line);
	@vcfRecord = split(/\t/,$line);
	$chrom = $vcfRecord[0];
	$vcfPos = $vcfRecord[1];
	$vcfRef = $vcfRecord[3];
	$vcfAlt = $vcfRecord[4];
	
	if ((length($vcfRef) == 1) && (length($vcfAlt) == 1)) {  #snp
		$start = $vcfPos;
		$end = $vcfPos;
		$ref = $vcfRef;
		$alt = $vcfAlt;		
	} elsif ( (length($vcfRef) > length($vcfAlt)) && (length($vcfAlt) == 1) )  {   #simple deletion
		$start = $vcfPos + length($vcfAlt);
		$end = $vcfPos + length($vcfRef) - length($vcfAlt);
		$ref = substr($vcfRef,length($vcfAlt));
		$alt = "-";
	} elsif ( (length($vcfRef) < length($vcfAlt)) && (length($vcfRef) == 1) ) {	#simple insertion
		$start = $vcfPos;
		$end = $vcfPos++;
		$ref = "-";
		$alt = substr($vcfAlt,length($vcfRef));
	} else {
	print "WARNING: vcf record is too complex for Oncotator web server: $line\n" and next;
	}

print TSV "$chrom\t$start\t$end\t$ref\t$alt\n";
}

close(VCF) or die "Unable to close input file: $!";
close(TSV) or die "Unable to close output file $!";


	


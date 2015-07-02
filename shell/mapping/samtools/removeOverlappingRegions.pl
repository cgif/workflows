#!/usr/bin/perl -w

use strict;

## script to generate interval list of non-overlapping amplicon regions
## input interval list need to be sorted by chromosome and position
## script does not check whether input is sorted

my $intList = $ARGV[0];

open(INPUT, $intList) or die "Unable to open input file: $!";
open (OUT, ">$ARGV[1]") or die "Unable to open output file: $!";

my $previous_line = "";

while (my $current_line = <INPUT>) {
	if ( $current_line =~ /^\@/ ) {
		print OUT $current_line;
	} elsif ($previous_line eq "") {
		$previous_line = $current_line;
		chomp $previous_line;
	} else {
		chomp $current_line;
		my @current_amplicon = split(/\t/, $current_line);
		my $current_chrom = $current_amplicon[0];
		my $current_start = $current_amplicon[1];
		my $current_end = $current_amplicon[2];
		my $current_strand = $current_amplicon[3];
		my $current_name = $current_amplicon[4];
		my @previous_amplicon = split(/\t/, $previous_line);
		my $previous_chrom = $previous_amplicon[0];
		my $previous_start = $previous_amplicon[1];
		my $previous_end = $previous_amplicon[2];
		my $previous_strand = $previous_amplicon[3];
		my $previous_name = $previous_amplicon[4];

		my $new_current_start;
		my $new_previous_end;

		if ( ( $current_chrom eq $previous_chrom )  && ( $current_start <= $previous_end ) ) {
			$new_current_start = ($previous_end+1);
			$new_previous_end = ($current_start-1);
			$previous_line = $previous_chrom . "\t" . $previous_start . "\t" . $new_previous_end . "\t" . "$previous_strand" . "\t" . $previous_name;
			$current_line = $current_chrom . "\t" . $new_current_start . "\t" . $current_end . "\t" . "$current_strand" . "\t" . $current_name;

			# print amplicon only if non-negative length (i.e. only if there is a region non-overlapping with any other amplicon)
#			if ( $previous_start < $new_previous_end ) {
#				print OUT $previous_line, "\n";
#			}

#		} else {
		}
		# print amplicon only if non-negative length (i.e. only if there is a region non-overlapping with any other amplicon)
		# checking here rather than in previous step, in case of amplicons within amplicons
		my @amplicon = split(/\t/, $previous_line);
		my $start = $amplicon[1];
		my $end = $amplicon[2];
		if ( $start < $end ) {		
			print OUT $previous_line, "\n";
		}

		$previous_line = $current_line;
	} 
}

print OUT $previous_line, "\n";



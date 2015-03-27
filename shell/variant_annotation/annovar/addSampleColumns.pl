#!/usr/bin/perl -w

use strict;

my $annovar = $ARGV[0];

if($annovar =~ /\.gz|\.gzip/){
	open(INPUT, "gzip -cd $annovar |") or die "Unable to open input file: $!";
} else {
	open(INPUT, $annovar) or die "Unable to open input file: $!";
}

open (OUT, ">$ARGV[1]") or die "Unable to open output file: $!";

my $header = <INPUT>;
print OUT "heterozygotes\thomozygotes\t$header";

$header =~ s/genotype_//g;

chomp $header;
my @header = split(/\t/,$header);

while (my $line = <INPUT>) {
	chomp $line;
	my @columns = split(/\t/,$line);
	my $search_het = "het";
	my( @indices_het )= grep { $columns[$_] eq $search_het } 0..$#columns;
	if (scalar @indices_het >= 1 ) { 
		foreach my $index_het ( @indices_het ) {
			print OUT "$header[$index_het],"
		}
	} else {
		print OUT "none";
	}
	print OUT "\t";
	my $search_hom = "hom";
	my( @indices_hom )= grep { $columns[$_] eq $search_hom } 0..$#columns;
	if (scalar @indices_hom >= 1 ) { 
		foreach my $index_hom ( @indices_hom ) {
			print OUT "$header[$index_hom],"
		}
	} else {
		print OUT "none";
	}
	print OUT "\t$line\n";
}

close (INPUT) or die "Unable to close input file: $!";
close (OUT) or die "Unable to close output file: $!";

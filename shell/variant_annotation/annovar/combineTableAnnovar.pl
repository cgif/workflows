#!/usr/bin/perl

use strict;


my $header;
my %pos2record;		# hash to contain variant annotations
my @genotypes;		# array of hashes containing all variants(keys) and genotypes and depths (values) for each sample.
my @samples;		# array to contain sample IDs

foreach my $file (@ARGV) {
	
	my @tokens = split(/\//,$file);
	my $file_name = @tokens[@tokens-1];
	my $sample = @tokens[@tokens-2];
	push @samples, $sample;
	my %sampleGenotype;		# hash which will go into @genotypes array
	
	if($file =~ /\.gz|\.gzip/){
		open(INPUT, "gzip -cd $file |") or die "Unable to open input file: $!";
	} else {
		open(INPUT, "<$file") or die "Unable to open input file: $!";
	}

	while(<INPUT>){
		chomp;
		my @columns = split(/\t/,$_);
		my $key = $columns[0]."_".$columns[1]."_".$columns[2]."_".$columns[3]."_".$columns[4];

		if($key eq "Chr_Start_End_Ref_Alt"){
			my @record = split(/\t/,$_);
			pop @record;				# remove "Other" field from avinput header		
			$header = join("\t", @record);
		} else {
			my $depth = pop @columns;
			my $qual = pop @columns;		
			my $genotype = pop @columns;
			push @columns, $qual;
			$pos2record{$key} = join("\t", @columns);
			$sampleGenotype{$key} = $genotype."\t".$depth;
		}
	} 	
	push @genotypes, {%sampleGenotype};
	close (INPUT) or die "Unable to close input file: $!";
}

#print header

#ANNOVAR header
#correct entries in ANNOVAR header
$header =~ s/\tbed\t/\tgenic_intolerance\t/;
$header =~ s/\tvcf\t/\tHGMD\t/;

print $header;
print "\tQUAL";

#sample genotype and depth header
foreach (@samples) {
	print "\tgenotype_".$_;
	print "\tdepth_".$_;
}
print "\n";
	
foreach my $key (sort (keys %pos2record)){			# by chromosome

	my $record = $pos2record{$key};				

#correct special entries for easier filtering 
#in genic intolerance scores
	$record =~ s/\tName=/\t/;
#in conservation and genomic duplications
	$record =~ s/\tScore=/\t/g;
	$record =~ s/;Name=\S+//g;

	print $record;	

	for (my $i = 0; $i < @ARGV; $i++) {
		if (exists $genotypes[$i]{$key}) {
		print "\t".$genotypes[$i]{$key};
		} else {
		print "\tNA\tNA";
		}
	}
print "\n";
}


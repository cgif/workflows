#!/usr/bin/perl
use strict;

my $header;
#my %pos2record;		# hash to contain variant annotations
#my @genotypes;		# array of hashes containing all variants(keys) and genotypes and depths (values) for each sample.
#my @samples;		# array to contain sample IDs

#my %exon2record;	# hash to contain exon annotations from cnv calls file (Ensembl gene, gene name)
my %cnvCalls;		# hash of hashes containing all called exons (keys) and calls (gain, loss, normal) for each sample.
my @samples;		# array to contain sample IDs
my %seen_samples;	
	
open(INPUT, $ARGV[0]) or die "Unable to open input file: $!";

while(<INPUT>){
	chomp;
	next if /^sample/;				#skip the header
	next if /^(\s)*$/;				#skip empty lines if present

	my @columns = split(/\t/,$_);
	my $sample = $columns[0];

	push(@samples, $sample) unless $seen_samples{$sample}++;

	

	my @exons = split(/,/,$columns[16]);
	
	foreach my $exon (@exons){
		my $exonKey = $exon;
		my $ratio = $columns[11];
		my $ensemblGene = $columns[15];
		my $geneName = $columns[17];
#		$exon2record{$exonKey} = $ensemblGene."\t".$geneName;
		$cnvCalls{$exonKey}{$sample} = $ratio;
	}
} 	

#get gene and position information
open(INGENE, $ARGV[1]) or die "Unable to open input file: $!";
my %exon2gene;

while(<INGENE>){
	chomp;
	next if /^Chromosome/;
	next if /^(\s)*$/;				#skip empty lines if present

	my @columns = split(/\t/,$_);
	my $exon = $columns[3];
	my $chr = $columns[0];
	my $start = $columns[1];
	my $end = $columns[2];
	my $ensgene = $columns[6];
	my $geneName = $columns[4];
	my $biotype = $columns[5];
	$exon2gene{$exon} = $ensgene."\t".$geneName."\t".$chr."\t".$start."\t".$end."\t".$biotype
} 	


#print header

print "exon\tGene\tGene_name\tChr\tStart\tEnd\tBiotype";

foreach my $sample (@samples){
	print "\t", $sample;
}
print "\n";

foreach my $key (keys %cnvCalls){
	if ($key ne "NA") {
		print $key;
		print "\t", $exon2gene{$key};

		foreach my $sample (@samples) {
			if (exists $cnvCalls{$key}{$sample}) {
				print "\t", $cnvCalls{$key}{$sample};
			} else {
			print "\tN";
			}
		}
	}
	print "\n";
}

close (INPUT) or die "Unable to close input file: $!";

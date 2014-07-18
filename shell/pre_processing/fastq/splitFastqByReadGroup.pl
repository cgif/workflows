#!/usr/bin/perl

#
# Splits an Illumina fastq file into lane read groups. 
# The input fastq file has to be read name sorted.
# Read group tags can be provided as input arguments 
# following the Illumina convention:
#
# <instrument>:<run number>:<flowcell ID>:<lane>
#
# If no read group tags provided they are parsed from the
# input file. 
#

use strict;

my $fastq_input = $ARGV[0];
#my @read_groups = ('HWI-ST759:244:C3R7BACXX:1','HWI-ST759:244:C3R7BACXX:2','HWI-ST759:244:C3R7BACXX:3','HWI-ST759:244:C3R7BACXX:4','HWI-ST759:244:C3R7BACXX:5','HWI-ST759:244:C3R7BACXX:6','HWI-ST759:244:C3R7BACXX:7');
my %file_handles;

print "splitting fastq file: $fastq_input\n";

#get read groups
my @read_groups;
for(my $i = 1; $i < @ARGV; $i++){
	push(@read_groups, $ARGV[$i]);
}

#if read groups are not provided by user
#parse them from file
if(@read_groups < 1){

	print "no read groups provided. Parsing from fastq file...\n";
	
	open(RG, "zcat $fastq_input | awk 'NR == 1 || (NR-1) % 4 == 0' | cut -f1,2,3,4 -d ':' | uniq | perl -pe 's/^@//' |") || die "Unable to open input fastq $fastq_input: $!";
	while(<RG>){
		push(@read_groups, $_);
	}
	close(RG);

	print "read groups: ";
	

} else {
	
	print "read groups provided: ";
	
}

foreach my $group (@read_groups) {
	print "$group ";
}		
print "\n";


#instanciate output file handles
foreach my $group (@read_groups) {
	$file_handles{$group} = anon_fh();
} 

#open output file handles
print "opening output file handles...";
foreach my $group (@read_groups) {
	my $fh = $file_handles{$group};
	my $filename = $group;
	$filename =~ s/:/_/g;
	$filename = $fastq_input."_".$filename.".fastq.gz";
	open($fh, "| gzip > $filename");
} 
print "done\n";

my $row_count = 0;
my $read_count = 0;

open (IN, "gzip -cd $fastq_input |") or die "Unable to open input fastq $fastq_input: $!";

while(<IN>){

	$read_count++;
	
	foreach my $group (@read_groups) {
		
		if(/$group/){
	
			my $fh = $file_handles{$group};
			#print header
			print { $fh } $_;
			#print sequence
			$_=<IN>;
			print { $fh } $_;
			#print header
			$_=<IN>;
			print { $fh } $_;
			#print quality
			$_=<IN>;
			print { $fh } $_;
			
			#break out of for loop
			last;

		} 	

	}

	if($read_count % 100000 == 0){
				print "$read_count reads processed\n"; 
	}	
	
}

close (IN);

#close output file handles
foreach my $group (@read_groups) {
	my $fh = $file_handles{$group};
	close($fh);
} 

sub anon_fh {
	local *FH;
	return *FH;
}

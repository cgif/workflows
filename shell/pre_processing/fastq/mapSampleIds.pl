#!/usr/bin/perl

#
# Renames sample directories in input directory from
# customer sample ID to facility sample ID

use strict;
use warnings;

my $input_directory = $ARGV[0];
my $read_group_info_file = $ARGV[1];
my $sample_id_mapping_file = $ARGV[2];


my %sample_facility_ids;
my %customer2facility_id;

#read facility sample IDs from read group
#info file
print "INFO: Reading read group info file $read_group_info_file\n";

open(RGI, $read_group_info_file) || die "ERROR: Unable to open $read_group_info_file: $!\n";

#skip header
<RGI>;

while(<RGI>){

	my @cols = split(/\t/);
	my $facility_id = $cols[1];
	if(defined $facility_id && $facility_id ne ""){
		$sample_facility_ids{$facility_id} = 0;
	}
}

close(RGI);

#read customer-to-facility ID mapping
print "INFO: Reading facility-to-customer sample ID mapping file $sample_id_mapping_file\n";

open(SIM, $sample_id_mapping_file) || die "ERROR: Unable to open $sample_id_mapping_file: $!\n";

#skip header
while(<SIM>){

	my @cols = split(/\t/);
	my $facility_id = $cols[0];
	my $customer_id = $cols[1];
	
	if(exists $sample_facility_ids{$facility_id}){
		
		$sample_facility_ids{$facility_id} = 1;
		$customer2facility_id{$customer_id} = $facility_id;
	}

}

close(SIM);

foreach my $facility_id (keys(%sample_facility_ids)){

	my $present = $sample_facility_ids{$facility_id};
	if($present == 0){
		print "WARNING: no mapping found for sample ID '$facility_id'.\n"
	}

}

#foreach my $customer_id (keys(%customer2facility_id)){
#
#	my $facility_id = $customer2facility_id{$customer_id};
#	
#	print "$customer_id\t$facility_id\n";
#
#}

print "INFO: Renaming sample folders in $input_directory\n";

opendir(INDIR, $input_directory) || die "ERROR: Unable to read $input_directory: $!\n";

my @files = readdir(INDIR);
	
foreach my $file (@files){

#print "$file\n";

	my $directory = "$input_directory/$file";

	if(-d $directory && !($file eq "." || $file eq "..")){

		my $customer_id = $file;

		#get facility ID
		if (exists $customer2facility_id{$customer_id}){
			
			my $facility_id = $customer2facility_id{$file};
			my $directory_renamed = "$input_directory/$facility_id";
				
			#rename directory
			my $command = "mv $directory $directory_renamed";
#			print "$command\n";
			my $exit_status = system($command); 
			if($exit_status != 0){
				die "ERROR: Failed to rename folder $directory.\n";
			}
		
		} else {
			print "WARNING: no facility ID found for sample $customer_id.\n"
		}
		
	}

}


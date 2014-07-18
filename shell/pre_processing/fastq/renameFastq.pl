#!/usr/bin/perl

#
# Renames fastq files in input directory to the format
# <run_id>_<barcode_sequence>_<lane>_<read_of_pair> and puts
# them into sample sub-directories. The name of the input
# directory must be the run ID. The sample sheet filename must
# have the format <run_id>.SampleSheet.csv. Fastq files must 
# have the file extension .fq or .fastq.

use strict;
use warnings;

my $input_directory = $ARGV[0];
my $sample_sheet = $ARGV[1];
my $platform = $ARGV[2];

#parse run ID from input directory
$input_directory=~s/\/$//;
opendir(INDIR, $input_directory) || die "ERROR: Unable to read $input_directory: $!";
my @tokens = split(/\//, $input_directory);
my $run_id = @tokens[@tokens-1];

#parse run ID from sample sheet
@tokens = split(/\//, $sample_sheet);
my $sample_sheet_file_name = $tokens[@tokens-1];
@tokens = split(/\./, $sample_sheet_file_name);
my $run_id_sample_sheet = $tokens[0];
open(SS, $sample_sheet) || die "ERROR: Unable to open $sample_sheet: $!";

#make sure input directory and sample sheet run ID match
if($run_id ne $run_id_sample_sheet){
	die "ERROR: Run ID of input directory and sample sheet do not match: $run_id != $run_id_sample_sheet";
}

my %sample_barcodes;
my %sample_names;

#read sample sheet

if($platform eq "miseq"){

	while(<SS>){

		if(/Sample_ID/){
	
			while(<SS>){
		
				my @cols = split(/,/);
				my $sample_id = $cols[0];
				my $sample_name = $cols[1];
				$sample_name =~ s/ /_/g;
				my $barcode = $cols[5];
				$barcode =~ s/ //g;
			
				$sample_barcodes{$sample_id} = $barcode;
				$sample_names{$sample_id} = $sample_name;
		
				#print "$sample_id\t$sample_name\t$barcode\n";
		
			}
		
		}

	}

} elsif ($platform eq "hiseq") {

	while(<SS>){

		if(/SampleID/){
	
			while(<SS>){
		
				my @cols = split(/,/);
				my $sample_id = $cols[2];
				$sample_id =~ s/ /_/g;
				my $barcode = $cols[5];
				$barcode =~ s/ //g;
				my $sample_lane = $cols[1];
				
				$sample_barcodes{$sample_id} = $barcode;
				$sample_names{$sample_id} = $sample_id
		
				#print "$sample_id\t$sample_name\t$barcode\n";
		
			}
		
		}

	}

} else {

	die "ERROR: Unsupported platform: $platform.\n";

}

print "renaming and reorganising files in directory: $input_directory\n";
print "run ID: $run_id\n";
print "platform: $platform\n";

#rename and move fastq files

my %processed_samples;

if($platform eq "miseq"){

	my @files = readdir(INDIR);
	
	foreach my $file (@files){

		if($file=~m/\.fq|\.fastq/){
		
			my @file_name_tokens = split(/_/, $file);
			my $sample_name = $file_name_tokens[0];
			my $sample_id = $file_name_tokens[1];
			$sample_id =~ s/S//;

			#get sample barcode
			if(exists $sample_barcodes{$sample_id}){
						
				$processed_samples{$sample_id} = 1;
				my $name = $sample_names{$sample_id};
				my $barcode = $sample_barcodes{$sample_id};
					
				my $sample_sub_dir = "$input_directory/$name";
				my $new_file_name = $run_id."_".$barcode;
				for(my $i = 2; $i < @file_name_tokens; $i++){
					$new_file_name = $new_file_name."_".$file_name_tokens[$i];
				}
				
				#create sample sub-folder
				my $command = "mkdir -p $sample_sub_dir";	 
				my $exit_status = system($command);
				if($exit_status != 0){
					die "ERROR: Failed to create sub-folder for sample $name: $command.\n";
				}
			
				#mv fastq file to sub-folder and rename
				$command = "mv $input_directory/$file ".$sample_sub_dir."/".$new_file_name;
				$exit_status = system($command); 
				if($exit_status != 0){
					die "ERROR: Failed to move fastq file $file to sub-folder: $command.\n";
				}
				
			} else {
				print "WARNING: No entry for sample $sample_name, ID $sample_id. Skipping file.";
			}		

		}
	
	}

} elsif ($platform eq "hiseq") {

	#read sample sub-directory names
	my @dirs = readdir(INDIR);

	#for each sample sub-directory
	foreach my $dir (@dirs){
	
		#split directory name and get sample name
		#without "Sample_" prefix
		my @dir_name_tokens = split(/_/, $dir);
		my $sample_name = $dir_name_tokens[1];
	
		#if sample name was parseable
		if (defined $sample_name){
		
			#check if sample in sample sheet...
			if(exists $sample_barcodes{$sample_name}){
		
				#add sample to list of processed samples
				$processed_samples{$sample_name} = 1;
			
				#get sample barcode
				my $barcode = $sample_barcodes{$sample_name};

				#mv sample sub-folder to directory without "Sample_" prefix
				my $sample_sub_dir = "$input_directory/$sample_name";
				my $command = "mv $input_directory/$dir $sample_sub_dir";
				my $exit_status = system($command);
				if($exit_status != 0){
					die "ERROR: Failed to move sample sub-folder $dir: $command.\n";
				}
		
				#get fastq files in sample sub-directory
				opendir(SUBINDIR, $sample_sub_dir) || die "ERROR: Unable to read sample sub-folder $input_directory/$sample_name: $!";
				my @files = readdir(SUBINDIR);
			
				#rename file by replacing the sample name with
				#the run ID
				foreach my $file (@files){
		
					if($file=~m/\.fq|\.fastq/ && $file=~m/$sample_name/ && $file=~m/$barcode/){
				
						my $destination_file = $file;
						$destination_file =~ s/$sample_name/$run_id/;
						
						my $command = "mv $sample_sub_dir/$file $sample_sub_dir/$destination_file";
						my $exit_status = system($command);
						if($exit_status != 0){
						die "ERROR: Failed to rename fastq file from $file to $destination_file: $command.\n";
						}
				
					}
		
				}
			
				close(SUBINDIR);
		
			#...skip sample if not in sample sheet and
			#throw warning
			} else {
					print "WARNING: No entry for sample $sample_name. Skipping files.";
			} # end of if(exists $sample_barcodes{$sample_name})
		
		} #end of if (defined $sample_name)
		
	} #end of foreach my $dir (@dirs)

} #end of if($platform eq "miseq")

close(INDIR);
close(SS);

foreach my $id (keys %sample_names){

	my $name = $sample_names{$id};
	if(!exists $processed_samples{$id}){
		print "WARNING: No fastq files found for record $name, ID $id in sample sheet.\n"; 
	}

}

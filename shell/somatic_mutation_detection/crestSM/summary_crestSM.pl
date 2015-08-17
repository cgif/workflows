#!/usr/bin/perl -w

use strict;
use warnings;

my $crest_dir="#inputDir";
my $deployment_server="#deploymentServer";
my $summary_deployment="#summaryDeployment";
my $project="#project";

#$crest_dir="/groupvol/cgi/results/johnson_glioma/crestSM/2015-01-19";
#$deployment_server="eliot.med.ic.ac.uk";
#$summary_deployment="/www/html/report/project/johnson_glioma/crestSM/2015-01-19";
#$project="johnson_glioma";

my $summary_link = $1 if  $summary_deployment =~ /\/www\/html\/(.*)/;
my(@samples, @data, $fraction, $start, $end, $start1, $end1, $start2, $end2, $color);

my $bed = "$crest_dir/multisample/$project.bed";
open(BED, ">$bed");
print BED "browser position chr1\n";
print BED "browser hide all\n";

opendir (INPUT, "$crest_dir") || print "Can't open $crest_dir\n";
while (defined(my $sample = readdir(INPUT))){

	next unless -e "$crest_dir/$sample/$sample.predSV.txt";
	push(@samples, $sample);
	print BED "track name=$sample type=bed visibility=3 db=hg19 itemRgb=On useScore=1\n";
	open (FILE, "$crest_dir/$sample/$sample.predSV.txt");

  	while (<FILE>){
		
		next if /hs37d5/;
		@data=split(/\t/,$_);  

		if ($data[8] =~ /INS|DEL|ITX|INV/){

			if ($data[1] <= $data[5]) {
				$start = $data[1]; 
				$end = $data[5];
			} else {
				$start = $data[5]; 
				$end = $data[1];
			} 
			
			$color = "255,0,0" if $data[8] eq "DEL"; 
			$color = "0,255,0" if $data[8] eq "INS"; 
			$color = "0,0,255" if $data[8] eq "ITX"; 
			$color = "0,0,255" if $data[8] eq "INV"; 
			
			print BED "chr$data[0]\t$start\t$end\t$data[8]\t0\t$data[2]\t$start\t$end\t$color\n";

		} elsif ($data[8] =~ /CTX/) {

			$start1 = $data[1]; $end1 = $start1 + 100;
			$start2 = $data[5]; $end2 = $start2 + 100;

			print BED "chr$data[0]\t$start1\t$end1\tCTX:chr$data[4]:$start2\t0\t$data[2]\t$start1\t$end1\t255,255,255\n";
			print BED "chr$data[4]\t$start2\t$end2\tCTX:chr$data[0]:$start1\t0\t$data[6]\t$start2\t$end2\t255,255,255\n";

 		}
	}
}

open (INDEX, ">$crest_dir/multisample/index.html");
print INDEX "<HTML><BR>\n";
print INDEX "<HEAD></HEAD><BR>\n";
print INDEX "<BODY><FONT SIZE = '+1'>Chromosomal rearangements discovered by CREST</FONT><BR>\n";
foreach my $sample (@samples){
	print INDEX "<A HREF = 'http://$deployment_server/$summary_link/$sample.predSV.html'>$sample</A><BR>\n";
}
print INDEX "<P><FONT SIZE = '+1'>Somatic chromosomal rearangements in <A HREF = 'http://$deployment_server/$summary_link/$project.bed'>BED</A> format<BR>\n"; 

system("scp -r $crest_dir/multisample/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("scp -r $crest_dir/multisample/$project.bed $deployment_server:$summary_deployment/$project.bed > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $summary_deployment/* > /dev/null 2>&1");


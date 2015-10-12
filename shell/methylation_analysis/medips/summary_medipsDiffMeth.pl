#!/usr/bin/perl -w

use strict;
use warnings;

my $medips_results_dir="#medipsResultsDir";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";
my $enc_dir = "#encDir";
my $random15 = "#random15";
my $sample_info = "#sampleInfo";

my $custom_tracks_URL = "$deployment_server/report/data/$random15/custom_track.txt";
open(TRACKS, ">$medips_results_dir/custom_track.txt");
print TRACKS "browser position chr21";

open(INDEX, ">$medips_results_dir/index.html");

my %data = ();
open(INFO, ">$sample_info");
while(<INFO>){

	chomp();
	/(\S+)\t(\S+)/;
	$data{$2}++;

}

foreach my $condition (keys %data){

	$bigWig_URL="$deployment_server/report/data/$random15/$condition.bw";
	print TRACKS "track type=bigWig name=$condition description=$condition.rpkm visibility=full db=hg19 viewLimits=0:10 autoScale=off windowingFunction=mean bigDataUrl=http://$bigWig_URL";

}

$URL="http://$deployment_server/ucsc/cgi-bin/hgTracks?org=human&db=hg19&position=chr21&hgct_customText=http://$custom_tracks_URL";
print INDEX "Methylation profiles can be seen at <A HREF=$URL>UCSC browser</A>";

system("scp $medips_results_dir/custom_track.txt $deployment_server:$enc_dir");
system("scp $medips_results_dir/index.html $deployment_server:$summary_deployment");
system("ssh $deployment_server chmod 0775 $summary_deployment/*");

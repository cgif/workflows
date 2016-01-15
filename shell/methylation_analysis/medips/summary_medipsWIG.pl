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
open(TRACKS, ">$medips_results_dir/WIG/custom_track.txt");
print TRACKS "browser position chr21\n";

open(INDEX, ">$medips_results_dir/WIG/index.html");
print INDEX "<TABLE CELLPADDING=5>";

my %data = ();
my $line = 1;
my ($bigWig_URL, $URL);
open(INFO, "$sample_info");
while(<INFO>){

	chomp();
	/(\S+)\t(\S+)/;
	if ($line == 1){
		print INDEX "<TR><TD><B>$1<TD><B>$2";
		$line++;
	}else{
		$data{$2}++;
		print INDEX "<TR><TD>$1<TD>$2";
	}

}
print INDEX "</TABLE>";

foreach my $condition (keys %data){

	next if $condition eq "category";
	$bigWig_URL="$deployment_server/report/data/$random15/$condition.bw";
	print TRACKS "track type=bigWig name=$condition description=$condition.rpkm visibility=full db=hg19 viewLimits=0:10 autoScale=off windowingFunction=mean bigDataUrl=http://$bigWig_URL\n";

}

$URL="http://$deployment_server/ucsc/cgi-bin/hgTracks?org=human&db=hg19&position=chr21&hgct_customText=http://$custom_tracks_URL";
print INDEX "<HR><P>Methylation profiles can be opened as custom tracks at <A HREF=$URL>UCSC browser</A>";

system("scp $medips_results_dir/WIG/custom_track.txt $deployment_server:$enc_dir");
system("scp $medips_results_dir/WIG/index.html $deployment_server:$summary_deployment");
system("ssh $deployment_server chmod 0775 $summary_deployment/*");
system("ssh $deployment_server chmod 0775 $enc_dir/*");

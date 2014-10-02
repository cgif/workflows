#!/usr/bin/perl -w
use strict;

#This script generates on-line summary that reflects jobs progress. 
#Once job is finished and non-empty log file is present in the run directory,
#If expected output file is not empty then green tick appears in the summary table.
#Otherwise red cross appears in the summary table.

my $project_dir_analysis = "projectDirAnalysis";
my $project_dir_results = "projectDirResults";
my $project = "#project";
my $date = "#today";
my $summary_results = "summaryResults";
my $deployment_server = "deploymentServer";
my $summary_deployment = "summaryDeployment";
my $data_deployment = "dataDeployment";
my $sample_list = "sampleList";
my $enc_dir = "encryptedDir";
my $summary_link = $1 if  $summary_deployment =~ /\/www\/html\/(.*)/;
my $data_link = $1 if $data_deployment =~ /\/www\/html\/(.*)/;


my @sample_all = ();

open (LIST, "$sample_list");
while (<LIST>){
	chomp();
	next if /^#/;
	next if /^(\s)*$/;
	my @record = split(/\t/,$_);
	push(@sample_all, $record[1]);
}

my @sample = Uniq(@sample_all);

my %sum = ();

foreach my $sample (@sample){
	my $callCNVs_log = "$project_dir_analysis/$date/$sample/run/${sample}.callCNVs.log";
	my $callCNVs_R_log= "$project_dir_analysis/$date/$sample/run/${sample}.callCNVs.R.log";
	if ((-s $callCNVs_log) && (-s $callCNVs_R_log)) {
	    my $R_grep = "";
	    $R_grep = `grep 'Execution halted' $callCNVs_R_log`;
		if ( $R_grep eq "" ){
			$sum{$sample}{'calling'} = "PASS";
	    }else{
			$sum{$sample}{'calling'} = "FAIL";
	    }

	    my $autosomes_cnvs = "$project_dir_results/$date/$sample/${sample}.cnv.calls.autosomes.tsv";
	    if (-s $autosomes_cnvs ){
			$sum{$sample}{'autosomes'} = "PASS";
	    }else{
			$sum{$sample}{'autosomes'} = "FAIL";
	    }

	    my $X_cnvs = "$project_dir_results/$date/$sample/${sample}.cnv.calls.X_chromosome.tsv";
	    if (-s $X_cnvs ){
			$sum{$sample}{'X_chromosome'} = "PASS";
	    }else{
			$sum{$sample}{'X_chromosome'} = "FAIL";
	    }
	}
}


open (OUT, ">$summary_results/index.html");

print OUT "<HTML>";
print OUT "<HEAD><META HTTP-EQUIV='refresh' CONTENT='60'></HEAD>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>CNV calling<BR>complete";
print OUT "<TH><CENTER>CNV calls<BR>autosomes";
print OUT "<TH><CENTER>CNV calls<BR>X chromosome";

#create data directory
system("ssh $deployment_server mkdir -p -m 775 $data_deployment/$enc_dir > /dev/null 2>&1");

my @tags = qw(autosomes X_chromosome);
foreach my $sample (sort @sample){
	print OUT "<TR><TD>$sample\n";

	if (defined $sum{$sample}{'calling'}){
		print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{'calling'} eq "PASS";
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{'calling'} eq "FAIL";
		print OUT "<TD> \n" if $sum{$sample}{'calling'} ne "FAIL" && $sum{$sample}{'calling'} ne "PASS";
	}
     
	foreach my $tag (@tags){
		if (defined $sum{$sample}{$tag}){
			my $cnv_calls = "$project_dir_results/$date/$sample/${sample}.cnv.calls.$tag.tsv";
			my $cnv_deployment = "$data_deployment/$enc_dir/${sample}.cnv.calls.$tag.tsv";
			if (-s $cnv_calls){
				system("scp -r $cnv_calls $deployment_server:$cnv_deployment > /dev/null 2>&1");
			}
				print OUT "<TD><CENTER><A HREF = http://$deployment_server/$data_link/$enc_dir/${sample}.cnv.calls.$tag.tsv><IMG BORDER='5' SRC=tick.png ALT=PASS></A>\n" if $sum{$sample}{$tag} eq "PASS";
				print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL></A>\n" if $sum{$sample}{$tag} eq "FAIL";
				print OUT "<TD> \n" if $sum{$sample}{$tag} ne "FAIL" && $sum{$sample}{$tag} ne "PASS";
		}
	}
}

print OUT "</TABLE>\n";
#########################################
#multisample results and summary metrics
#########################################
#horizontal line
print OUT "<HR>";

#deploy files
my $multisample_calls = "$project_dir_results/$date/multisample/$project.cnvs.tsv";
my $multisample_bed = "$project_dir_results/$date/multisample/$project.cnvs.bed";
my $multisample_summary = "$project_dir_results/$date/multisample/$project.cnvs.summary.tsv";
my $multisample_summary_php = "$project_dir_results/$date/multisample/$project.cnvs.summary.php";
my $multisample_exons = "$project_dir_results/$date/multisample/$project.cnvs.exons.tsv";

#create target directory on server
system("ssh $deployment_server mkdir -p -m 775 $summary_deployment/multisample > /dev/null 2>&1");
system("ssh $deployment_server mkdir -p -m 775 $data_deployment/$enc_dir/multisample > /dev/null 2>&1");




#copy files to server
system("scp -r $multisample_calls $deployment_server:$data_deployment/$enc_dir/multisample > /dev/null 2>&1");
system("scp -r $multisample_bed $deployment_server:$data_deployment/$enc_dir/multisample > /dev/null 2>&1");
system("scp -r $multisample_summary $deployment_server:$summary_deployment/multisample > /dev/null 2>&1");
system("scp -r $multisample_summary_php $deployment_server:$summary_deployment/multisample > /dev/null 2>&1");
system("scp -r $multisample_exons $deployment_server:$data_deployment/$enc_dir/multisample > /dev/null 2>&1");

my $multisample_calls_url = "http://$deployment_server/$data_link/$enc_dir/multisample/$project.cnvs.tsv";
my $multisample_bed_url = "http://$deployment_server/$data_link/$enc_dir/multisample/$project.cnvs.bed";
my $multisample_summary_url = "http://$deployment_server/$summary_link/multisample/$project.cnvs.summary.tsv";
my $multisample_summary_php_url = "http://$deployment_server/$summary_link/multisample/$project.cnvs.summary.php";
my $multisample_exons_url = "http://$deployment_server/$data_link/$enc_dir/multisample/$project.cnvs.exons.tsv";

print OUT "<P><FONT SIZE = '+1'><A HREF = '$multisample_calls_url'>Multisample CNV calls [TSV]</A> ".
							   "<A HREF = '$multisample_bed_url'>[BED]</A> ".
							   "<A HREF = '$multisample_exons_url'>[EXONS]</A> ".
			  "</FONT><BR>";
			  
print OUT "<P><FONT SIZE = '+1'><A HREF = '$multisample_summary_php_url'>View summary statistics</A> ".
							   "<A HREF = '$multisample_summary_url'>[TSV]</A> ".
			"</FONT><BR>";


####################
#build UCSC track - doesn't make sense, because we don't have CNV cb tracks
####################

#system("scp -r $bed $deployment_server:$summary_deployment/multisample > /dev/null 2>&1");

#my $customTrackURL = "http://$deployment_server/$summary_link/multisample/exomeDepth_test.cnvs.bed";
#my $directoryURL = "http://$deployment_server/$summary_link/multisample/";
#my $customTrack = "track type=bigBed name=$project visibility=full bigDataUrl=$customTrackURL";
#my $url = "http://$deployment_server/ucsc/cgi-bin/hgTracks?org=human&db=hg19&position=chr21:33,031,597-33,041,570&hgct_customText=$customTrack";

#if (defined $sum{'multisample'}{'raw_vcf'}){
#    if ($sum{'multisample'}{'raw_vcf'} eq "PASS"){
#	print OUT "<P><FONT SIZE = '+1'>Results can be viewed as a <A HREF = '$url'>custom track</A> in the UCSC Genome Browser</FONT><BR>";
#	print OUT "<P><FONT SIZE = '+1'>Download vcf file <A HREF = '$directoryURL'>here</A></FONT><BR>";
#    }
#}


#####################
#copy index.html
#####################

system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $summary_deployment/* > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 775 $data_deployment/$enc_dir/* > /dev/null 2>&1");

close (OUT);

## Subroutine to get unique sample names

sub Uniq{
	my %temp_hash = map { $_, 0 } @_;
	return keys %temp_hash;
}


#!/usr/bin/perl -w

use strict;
use warnings;

my $fastq_dir = "#pathReadsFastq";
my $analysis_dir = "#pathAnalysisDir";
my $report_dir = "#pathReportsDir";
my $ms_report_dir = "#pathMSReportsDir";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

#collect data from log files on mapping statistics
my %sum = ();

#open rawdata project directory
opendir (PROJECT, "$report_dir") || print "Can't open $report_dir\n";

#for each sample in the fastq directory... 
while (defined(my $sample = readdir(PROJECT))){

	#skip . and .. and multisample directory
    next if ($sample =~ /^\./ || $sample =~ /^multisample$/);
    
    #assert that path is a directory path
    #skip if not
    my $sample_dir = "$report_dir/$sample";
    next unless (-d "$sample_dir");

    #check for miRDeep run log file
    my $run_dir = "$analysis_dir/$sample/run";
    opendir (RUNDIR, "$run_dir");
    
    my $log_count = 0;
    while (defined(my $log = readdir(RUNDIR))) {
    	
    	#skip . and ..
    	next if $log =~ /^\./ || !($log =~ /\.log/);
    
    	#count log files
    	$log_count++;
    	
    	#parse errors
    	$sum{$sample}{'error'}=`grep ERROR $run_dir/$log`;

		#get mapping statistics
		$sum{$sample}{'total'}=`grep 'seq:' $run_dir/$log | tail -n +2 | cut -f 1 | cut -d ' ' -f 2`;
		$sum{$sample}{'mapped'}=`grep 'seq:' $run_dir/$log | tail -n +2 | cut -f 2`;
		$sum{$sample}{'unmapped'}=`grep 'seq:' $run_dir/$log | tail -n +2 | cut -f 3`;
		$sum{$sample}{'prMapped'}=`grep 'seq:' $run_dir/$log | tail -n +2 | cut -f 4`;
		$sum{$sample}{'prUnmapped'}=`grep 'seq:' $run_dir/$log | tail -n +2 | cut -f 5`;

	}
	
	#if no run logs existed for sample... 
	if($log_count == 0){
		$sum{$sample}{'error'}="miRDeep2 analysis was not done for this sample.";
	}
	
	#skip if errors occured
    next if $sum{$sample}{'error'} =~ /\w+/;
    
}


#print extracted data into summary file

#open output file
open (OUT, ">$ms_report_dir/index.html");

print OUT "<HTML>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Expressed miRNAs";
print OUT "<TH><CENTER>Results file";
print OUT "<TH><CENTER>Total reads";
print OUT "<TH><CENTER>Mapped reads";
print OUT "<TH><CENTER>Unmapped reads";
print OUT "<TH><CENTER>Proportion<BR>mapped";
print OUT "<TH><CENTER>Proportion<BR>unmapped\n";

my @tags = qw(total mapped unmapped prMapped prUnmapped);

#for each sample (alphanumerically sorted)
foreach my $sample (sort {$a cmp $b} keys %sum){

    print OUT "<TR><TD>$sample";
    
	if ($sum{$sample}{'error'} =~ /\w+/){

		print OUT "<TD COLSPAN=10>$sum{$sample}{'error'}";

	} else {
	    
		my $link_html = $summary_deployment;
		$link_html =~ s/www\/html\/(.*)$/$1/;
		$link_html = "http://"."$deployment_server/"."$link_html/"."$sample/$sample.expression.html";
		
		print OUT "<TD><CENTER><A HREF = $link_html>Expression profile</A>";

		my $link_tsv = $summary_deployment;
		$link_tsv =~ s/www\/html\/(.*)$/$1/;
		$link_tsv = "http://"."$deployment_server/"."$link_tsv/"."$sample/$sample.miRNA_expressed.tsv";
		
		print OUT "<TD><CENTER><A HREF = $link_tsv>Expression table [TSV]</A>";

		foreach my $tag (@tags){
			print OUT "<TD><CENTER>$sum{$sample}{$tag}";
		}
	}
		
	print OUT "\n";

}

print OUT "</TABLE></HTML>";

system("chmod 660 $ms_report_dir/index.html");
system("scp -r $ms_report_dir/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");



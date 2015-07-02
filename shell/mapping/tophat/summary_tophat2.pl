#!/usr/bin/perl -w

use strict;
use warnings;

my $project_dir_analysis = "#projectDirAnalysis";
my $project_dir_results = "#projectDirResults";
my $summary_results = "#summaryResults";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

$project_dir_analysis =~ s/^(\S+)\/\S+$/$1/;
$project_dir_results =~ s/^(\S+)\/\S+$/$1/;

#collect data from bwa log files
my %sum = (); 

#collect info for each sample run
my($sample, $sample_dir, $script, $output_prefix, $read, $bam, $log, $summary, $grep_mapping, $grep_concordant, $left, $right);

opendir (PROJECT, "$project_dir_analysis"); 
while (defined($sample = readdir(PROJECT))){
    next if $sample =~ /^\./;

    $sample_dir = "$project_dir_analysis/$sample/run";
    next unless (-d "$sample_dir");

    #collect info for each read group
    opendir (SAMPLE, "$sample_dir");
    while (defined($script = readdir(SAMPLE))){
        next unless $script =~ /tophat2.((.*)\.vs\..*).sh/;
	$output_prefix = $1;
	$read = $2;

        #check that the tophat run is finished, bam is generated and log file has no warning messages
	$bam = "$project_dir_results/$sample/$output_prefix".".sorted.bam";
	$log = "$sample_dir/tophat2."."$output_prefix".".log";

        if (-e "$bam"){
            $sum{$sample}{$read}{'bam'} = "PASS";
        }

	if (-e "$log"){
            if (`grep 'Skipped' $log`){
                $sum{$sample}{$read}{'mate'} = "FAIL";
            }else{
                $sum{$sample}{$read}{'mate'} = "PASS";
            }

            if (`grep -i 'WARNING' $log|grep -v '20bp'|grep -v 'picard'`){
                $sum{$sample}{$read}{'bam'} = "FAIL";
	    }  
        }

	# collect statistics from tophat run
	$summary = "$project_dir_results/$sample/$output_prefix".".align_summary.txt";
	if (-e $summary){
            $grep_mapping = "";
	    $grep_mapping = `grep 'overall read mapping rate' $summary`;
	    if ($grep_mapping =~ /mapping/){
                $grep_mapping =~ s/^(.*) overall read mapping rate./$1/;
		$sum{$sample}{$read}{'mapping'} = "$grep_mapping";
	    }

	    $grep_concordant = "";
	    $grep_concordant = `grep 'concordant pair alignment rate' $summary`;
	    if ($grep_concordant =~ /concordant/){
                $grep_concordant =~ s/^(.*) concordant pair alignment rate./$1/;
		$sum{$sample}{$read}{'concordant'} = "$grep_concordant";
	    }

	    $left = 0; $right = 0;
	    open (SUMMARY, "$summary");
	    while(<SUMMARY>){
	        $left++ if /Left reads/;
		$left = 0 if /Right reads/;
		$right++ if /Right reads/;
		$right = 0 if /Unpaired reads/;
		$sum{$sample}{$read}{'in'} = $1 if /Input\s+:\s+(\d+)/ && $left;
		$sum{$sample}{$read}{'mapped'} = $1 if /Mapped\s+:\s+(\d+)/ && $left;
		$sum{$sample}{$read}{'mult'} = $1 if /of these:\s+(\d+)/ && $left;
		$sum{$sample}{$read}{'in'} += $1 if /Input\s+:\s+(\d+)/ && $right;
		$sum{$sample}{$read}{'mapped'} += $1 if /Mapped\s+:\s+(\d+)/ && $right;
		$sum{$sample}{$read}{'mult'} += $1 if /of these:\s+(\d+)/ && $right;
            }
	}
    }
}

#print extracted data into summary file
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Fastq files";
print OUT "<TH><CENTER>Mate files";
print OUT "<TH><CENTER>Alignment";
print OUT "<TH><CENTER>Input reads";
print OUT "<TH><CENTER>Mapped reads";
print OUT "<TH><CENTER>Reads with multiple alignments";
print OUT "<TH><CENTER>Overall Read Mapping Rate";
print OUT "<TH><CENTER>Concordant Pair Alignment Rate\n";

foreach my $sample (sort {$a cmp $b} keys %sum){
    print OUT "<TR><TD>$sample";
    my $count = 0;

    foreach my $read (sort {$a cmp $b} keys %{$sum{$sample}}){
	print OUT "<TD>$read" if $count == 0;
	print OUT "<TR><TD><TD>$read" if $count;
	$count++;

	if (defined($sum{$sample}{$read}{'mate'})){
            if ($sum{$sample}{$read}{'mate'} eq "PASS"){
    	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
	    }elsif ($sum{$sample}{$read}{'mate'} eq "FAIL"){
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	    }
	}else{
	    print OUT "<TD><CENTER>";
	}

	if (defined($sum{$sample}{$read}{'bam'})){
	    if ($sum{$sample}{$read}{'bam'} eq "PASS"){
		print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
	    }elsif ($sum{$sample}{$read}{'bam'} eq "FAIL"){
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	    }
	}else{
	    print OUT "<TD><CENTER>";
	}

	if (defined($sum{$sample}{$read}{'in'})){
	    print OUT "<TD><CENTER>$sum{$sample}{$read}{'in'}";
	    print OUT "<TD><CENTER>$sum{$sample}{$read}{'mapped'}";
	    print OUT "<TD><CENTER>$sum{$sample}{$read}{'mult'}";
	    print OUT "<TD><CENTER>$sum{$sample}{$read}{'mapping'}";
	    print OUT "<TD><CENTER>$sum{$sample}{$read}{'concordant'}";
	}
    }
}
print OUT "</TABLE></HTML>";

#copy summary file on eliot 
system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");


#!/usr/bin/perl -w

use strict;
use warnings;

my $project_dir_analysis = "#projectDirAnalysis";
my $project_dir_results = "#projectDirResults";
my $summary_results = "#summaryResults";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

#my $project_dir_analysis = "/project/tgu/runs/sharma_methynets-amplicon/bismark/2016-12-04/multisample";
#my $project_dir_results = "/project/tgu/results/sharma_methynets-amplicon/bismark/2016-12-04/multisample";
#my $summary_results = "/project/tgu/runs/sharma_methynets-amplicon/bismark/2016-12-04/multisample";
#my $deployment_server = "eliot.med.ic.ac.uk";
#my $summary_deployment = "/www/html/report/project/sharma_methynets-amplicon/bismark/2016-12-04";


$project_dir_analysis =~ s/^(\S+)\/\S+$/$1/;
$project_dir_results =~ s/^(\S+)\/\S+$/$1/;

#collect data from bismark report files
my %sum = (); 

#collect info for each sample run
#my($sample, $sample_dir, $script, $output_prefix, $read, $bam, $log, $summary, $grep_mapping, $grep_concordant, $left, $right);

my ($sample, $sample_dir, $script, $bam, $log, $report, $grep_total, $grep_unique, $grep_eff, $grep_unmapped, $grep_non_un, $grep_discarded, $grep_top, $grep_bottom);

opendir (PROJECT, "$project_dir_analysis"); 
while (defined($sample = readdir(PROJECT))){
    next if $sample =~ /^\./;

    $sample_dir = "$project_dir_analysis/$sample/run";
    next unless (-d "$sample_dir");

    #collect info for each read group
    opendir (SAMPLE, "$sample_dir");
    while (defined($script = readdir(SAMPLE))){

        next unless $script =~ /BSIGF(.*)sh/;

        #check that the bismark run is finished and bam is generated and log file has no error messages
	$bam = "$project_dir_results/$sample/bam/${sample}_pe.bam";

	$log = "$sample_dir/BS"."$sample".".log";

        if (-e "$bam"){
            $sum{$sample}{'bam'} = "PASS";
        }

	if (-e "$log"){
            if (`grep 'ERROR' $log`){
                $sum{$sample}{'bam'} = "FAIL";
            }else{
                $sum{$sample}{'bam'} = "PASS";
            }
	}

	# collect statistics from bismark run
	my $summary = "$project_dir_results/$sample/bam/$sample"."_PE_report.txt";
	if (-e $summary){
            $grep_total = "";
	    $grep_total = `grep 'Sequence pairs analysed in total' $summary`;
	    if ($grep_total =~ /analysed/){
                $grep_total =~ s/.*\t(\d*)$/$1/;
		$sum{$sample}{'total'} = "$grep_total";
	    }

	    $grep_unique = "";
	    $grep_unique = `grep 'Number of paired-end alignments with a unique best hit' $summary`;
	    if ($grep_unique =~ /alignments/){
                $grep_unique =~ s/.*\t(\d*)$/$1/;
		$sum{$sample}{'unique'} = "$grep_unique";
	    }

	    $grep_eff = "";
	    $grep_eff = `grep 'Mapping efficiency' $summary`;
	    if ($grep_eff =~ /efficiency/){
                $grep_eff =~ s/.+\t(\d*)/$1/;
		$sum{$sample}{'eff'} = "$grep_eff";
	    }

	    $grep_unmapped = "";
	    $grep_unmapped = `grep 'Sequence pairs with no alignments under any condition' $summary`;
	    if ($grep_unmapped =~ /Sequence/){
                $grep_unmapped =~ s/.*\t(\d*)$/$1/;
		$sum{$sample}{'unmapped'} = "$grep_unmapped";
	    }

	    $grep_non_un = "";
	    $grep_non_un = `grep 'Sequence pairs did not map uniquely' $summary`;
	    if ($grep_non_un =~ /Sequence/){
                $grep_non_un =~ s/.*\t(\d*)$/$1/;
		$sum{$sample}{'non_unique'} = "$grep_non_un";
	    }

	    $grep_discarded = "";
	    $grep_discarded = `grep 'Sequence pairs which were discarded' $summary`;
	    if ($grep_discarded =~ /discarded/){
                $grep_discarded =~ s/.*\t(\d*)$/$1/;
		$sum{$sample}{'discarded'} = "$grep_discarded";
	    }

	    $grep_top = "";
	    $grep_top = `grep 'CT/GA/CT' $summary`;
	    if ($grep_top =~ /converted/){
                $grep_top =~ s/\D+(\d+)\t\(\(converted\) top strand\)$/$1/;
		$sum{$sample}{'top'} = "$grep_top";
	    }

	    $grep_bottom = "";
	    $grep_bottom = `grep 'CT/GA/GA' $summary`;
	    if ($grep_bottom =~ /converted/){
                $grep_bottom =~ s/\D+(\d+)\t\(\(converted\) bottom strand\)$/$1/;
		$sum{$sample}{'bottom'} = "$grep_bottom";
	    }

        } ## end of if (-e $summary)
    } ## end of  while (defined($script = readdir(SAMPLE)))
} ## end of while (defined($sample = readdir(PROJECT)))

#print extracted data into summary file
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>BAM files";
print OUT "<TH><CENTER>Total read pairs";
print OUT "<TH><CENTER>Uniquely aligned";
print OUT "<TH><CENTER>Mapping efficiency";
print OUT "<TH><CENTER>Unmapped pairs";
print OUT "<TH><CENTER>Pairs with multiple alignments";
print OUT "<TH><CENTER>Discarded pairs";
print OUT "<TH><CENTER>Aligned to top strand";
print OUT "<TH><CENTER>Aligned to bottom strand\n";

foreach my $sample (sort {$a cmp $b} keys %sum){
	print OUT "<TR><TD>$sample";
#    my $count = 0;

#    foreach my $read (sort {$a cmp $b} keys %{$sum{$sample}}){
#	print OUT "<TD>$read" if $count == 0;
#	print OUT "<TR><TD><TD>$read" if $count;
#	$count++;

	if (defined($sum{$sample}{'bam'})){
            if ($sum{$sample}{'bam'} eq "PASS"){
    	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
	    }elsif ($sum{$sample}{'bam'} eq "FAIL"){
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	    }
	}else{
	    print OUT "<TD><CENTER>";
	}

	if (defined($sum{$sample}{'total'})){
	    print OUT "<TD><CENTER>$sum{$sample}{'total'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'unique'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'eff'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'unmapped'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'non_unique'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'discarded'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'top'}";
	    print OUT "<TD><CENTER>$sum{$sample}{'bottom'}";
	}
}
print OUT "</TABLE></HTML>";

#copy summary file on eliot 
system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");


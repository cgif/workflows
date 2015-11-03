#!/usr/bin/perl -w

use strict;
use warnings;

my $medips_results_dir="#medipsResultsDir";
my $medips_runs_dir="#medipsRunsDir";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

#for each sample
my @samples = ();
my %data = ();
my ($log);

opendir(DIR, $medips_results_dir) or die $!;
    while (my $sample = readdir(DIR)) {

	next if $sample =~ /multisample|\./;
	push(@samples, $sample);
	$log = "$medips_runs_dir/$sample/$sample.qualityMetrics.log";
	open(LOG, "$log") || print "Can't open $log\n";
	while(<LOG>){
		$data{$sample}{'total'} = int($1/1000000) if /#estimated correlation for artifically doubled set of reads\s+(\d+) \d+/;
		$data{$sample}{'est_cor'} = $1 if /#estimated correlation for artifically doubled set of reads\s+\d+\s+(\d+\.\d\d)/;
		$data{$sample}{'no_GC'} = $1/($data{$sample}{'total'}*10000) if /#number of reads which do not cover CpG\s+(\d+)/;
		$data{$sample}{'no_GC'} =~ s/(\d+\.\d\d)\d+/$1/ if /#number of reads which do not cover CpG\s+(\d+)/;
		$data{$sample}{'CG_no_cov'} = $1 if /#number of CpG that are not covered\s+(\d+\.\d\d)/;
		$data{$sample}{'enrich_GC'} = $1 if /#observed\/expected ratio of CpGs within the covered regions \/ same ratio within the genome\s+(\d+\.\d\d)/; 

	}
	system("scp $medips_results_dir/$sample/$sample.sat_plot.png $deployment_server:$summary_deployment/$sample.sat_plot.png > /dev/null 2>&1");
	system("scp $medips_results_dir/$sample/$sample.seqcov_pie.png $deployment_server:$summary_deployment/$sample.seqcov_pie.png > /dev/null 2>&1");

}


#print summary
open (OUT, ">$medips_results_dir/multisample/index.html");
print OUT "<HTML>";
print OUT "<HEAD><TITLE>MEDIPS - SAMPLES QUALITY CONTROL</TITLE></HEAD>";
print OUT "<BODY>";
print OUT "<P>*Saturation analysis estimate, whether the number of short reads is sufficient to generate a saturated and reproducible coverage profile. The value is the Pearson correlation obtained by considering the artificially doubled set of reads (should be close to 1)<BR>\n";
print OUT "<P>**CpG enrichment value shows how strong the covered regions are enriched for CpGs compared to the reference genome (observed vs expected ratio of CpGs within the covered regions divided by same ratio for the entire genome (should be above 1, input samples have value in the range of 1.0 - 1.1)<BR>\n";
print OUT "<P>The better your enrichment experiment worked the less reads without CpGs you see (input samples usually show above 10%). Also, the fraction of uncovered genomic CpGs (x0) as well as the fraction of genomic CpGs with high coverage (>x5) should increase following the enrichment.<BR>\n";
print OUT "<P><HR><TABLE CELLPADDING=5><TR>";

my $sample_count = 0;
foreach my $sample (@samples){

	print OUT "<TR><TD COLSPAN = 5><HR></TABLE><TABLE CELLPADDING=5><TR>" if $sample_count % 5 == 0 && $sample_count > 0;
	print OUT "<TD><FONT SIZE = +1><CENTER><B>$sample</B></CENTER>\n";
	print OUT "<BR><CENTER>total read pairs $data{$sample}{'total'} M</CENTER>\n";
	print OUT "<P><IMG WIGTH=300 HEIGHT=300 SRC=$sample.sat_plot.png>\n";
	print OUT "<BR><CENTER>saturation* $data{$sample}{'est_cor'}</CENTER>\n";
	print OUT "<P><IMG WIGTH=300 HEIGHT=300 SRC=$sample.seqcov_pie.png>\n";
	print OUT "<BR><CENTER>no CpG reads $data{$sample}{'no_GC'} %</CENTER>\n";
	print OUT "<BR><CENTER>uncovered CpGs $data{$sample}{'CG_no_cov'} %</CENTER>\n";
	print OUT "<BR><CENTER>CpG enrichement** $data{$sample}{'enrich_GC'}</CENTER>\n";
	$sample_count ++;

}
print OUT "</TABLE>";

my $cor_matrix="$medips_results_dir/multisample/correlation_matrix.tsv";
if (-e $cor_matrix){

	system("echo -ne '\t'|cat - $medips_results_dir/multisample/correlation_matrix.tsv > $medips_results_dir/multisample/tmp.tsv");
	system("mv $medips_results_dir/multisample/tmp.tsv $medips_results_dir/multisample/correlation_matrix.tsv");
	print OUT "<HR><P><A HREF = correlation_matrix.tsv>Correlation Matrics</A> and <A HREF = correlation_plot.png>Hierarchical Clustering</A> for genome coverage profiles across all samples\n";
	system("scp $medips_results_dir/multisample/correlation_matrix.tsv $deployment_server:$summary_deployment/correlation_matrix.tsv > /dev/null 2>&1");
	system("scp $medips_results_dir/multisample/correlation_plot.png $deployment_server:$summary_deployment/correlation_plot.png > /dev/null 2>&1");

}

system("scp $medips_results_dir/multisample/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod 0664 $summary_deployment/* > /dev/null 2>&1");

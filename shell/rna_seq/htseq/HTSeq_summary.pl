#!/usr/bin/perl -w

use strict;
use warnings;

my $input_dir = "#inputDir";
my $project_dir_results = "#projectDirResults";
my $summary_results = "#summaryResults";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

#collect data from htseq output files for each sample
my %sum = (); 

opendir (PROJECT, "$input_dir"); 
while (defined(my $sample = readdir(PROJECT))){

    next if $sample =~ /^\./;
    next if $sample eq "multisample";

    my $flagstat = "$input_dir/$sample/$sample.bam.flagstat";
    if (-f "$flagstat"){
        $sum{$sample}{'total'} = `grep total $flagstat`;
        $sum{$sample}{'total'} =~ s/^(\d+)\s.*/$1/;
    }

    my $table = "$project_dir_results/$sample/$sample.HTSeq.counts.discard";
    next unless (-f "$table");
    foreach my $param qw(no_feature ambiguous too_low_aQual not_aligned alignment_not_unique){
        $sum{$sample}{$param} = `tail $table|grep '$param'`;
        $sum{$sample}{$param} =~ s/\w+\s+(\d+)/$1/;
    }
}


#print extracted data into summary file
open (OUT, ">$summary_results/index.html");
print OUT "<HTML>";
print OUT "<TABLE CELLPADDING=5><TR>";
print OUT "<TH ROWSPAN = 2><CENTER>Sample";
print OUT "<TH ROWSPAN = 2><CENTER>Total reads";
print OUT "<TH COLSPAN = 5><CENTER>Discarded reads";
print OUT "<TR><TH><CENTER>Reads that don't overlap with any features";
print OUT "<TH><CENTER>Reads that overlap with > 1 feature";
print OUT "<TH><CENTER>Reads with alignment quality below 10";
print OUT "<TH><CENTER>Reads that did not align";
print OUT "<TH><CENTER>Reads that align to multiple locations";
print OUT "\n";

foreach my $sample (keys %sum){

    print OUT "<TR><TD>$sample"; 
    print OUT "<TD><CENTER>$sum{$sample}{'total'}"; 

    foreach my $param qw(no_feature ambiguous too_low_aQual not_aligned alignment_not_unique){
        print OUT "<TD><CENTER>$sum{$sample}{$param}";
	if (defined($sum{$sample}{'total'})){
	    my $percent = ($sum{$sample}{$param}/$sum{$sample}{'total'})*100;
	    $percent =~ s/(\d+\.\d\d).*/$1/;
            print OUT "($percent%)";
	}
    }

    print OUT "\n";

}
system("chmod 0660 $summary_results/index.html");

#copy summary file on eliot 
system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");

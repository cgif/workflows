#!/usr/bin/perl -w

use strict;
use warnings;

my $fastq_dir = "#pathReadsFastq";
my $analysis_dir = "#pathAnalysisDir";
my $report_dir = "#pathReportsDir";
my $ms_report_dir = "#pathMSReportsDir";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

#collect data from fastqc summaries in %sum
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

    #check for FastQC run log file; integrity and proper pairing of fastq files
    my $run_dir = "$analysis_dir/$sample";
    opendir (RUNDIR, "$run_dir");
    
    my $log_count = 0;
    while (defined(my $log = readdir(RUNDIR))) {
    	
    	#skip . and ..
    	next if $log =~ /^\./ || !($log =~ /\.log/);
    
    	#count log files
    	$log_count++;
    	
    	#parse errors
    	$sum{$sample}{'error'}=`grep ERROR $run_dir/$log`;

    }
	
    #if no run logs existed for sample... 
    if($log_count == 0){

	$sum{$sample}{'error'}="FastQC check was not done for this sample.";

    }
	
    #skip if errors occured
    next if $sum{$sample}{'error'} =~ /\w+/;

    #for each fastq retrieve FAIL/WARN/PASS for each quality metrics
    
    #open sample fastq directory
    opendir (SAMPLE, "$sample_dir") || print "Can't open $sample_dir\n";
    
    #for each fastqc report
    while (defined(my $fastqc_report = readdir(SAMPLE))){
    
    	#check fastqc extension
	next unless $fastqc_report  =~ /^(\S+)_fastqc$/;
			
	#if fastqc summary report exists...
        if (-e "$report_dir/$sample/$fastqc_report/summary.txt"){
	    	
	    my $cat = 11;
	    	
	    #...parse fastqc summary report
	    open (SUM, "$report_dir/$sample/$fastqc_report/summary.txt") || print "Can't open $report_dir/$sample/$fastqc_report/summary.txt\n";
	    while (<SUM>){
        
		next if /Basic Statistics/;
	    
	    	if (/^(\w\w\w\w)\t/){

		    $sum{$sample}{$fastqc_report}{$cat} = $1;
	    	    $cat++;

		}
				
	    }
			
	} else {
		
            $sum{$sample}{$fastqc_report}{'error'}="FastQC check was not done for this fastq file.";
	    next;
	    	
	}

        #read positions for low quality bases and total sequences produced for each read for each sample
	my $f1 = 0; 
	my $f3 = 0;
        my @l_qual = ();
        my $l_qual = ""; 
        my $seq_tot = 0;
	
	#open fastqc data file
	open (DATA, "$report_dir/$sample/$fastqc_report/fastqc_data.txt")||print "Can't open $report_dir/$sample/$fastqc_report/fastqc_data.txt\n";
		
	while (<DATA>){
	    	
	    $f1 = 1 if /^>>Per base sequence quality/;
	   
	    if ($f1 && /^(\d+)\t\S+\t(\d+)\S+\t(\d+)\S+\t/){

            	my $pos = $1;
           	my $med = $2;
            	my $lq = $3;
            	push(@l_qual, $pos) if ($med < 20 || $lq < 5);     
               
	    }
	    	
	    $f1 = 0 if /^>>END_MODULE/;
	    $f3 = 1 if /^>>Basic Statistics/;
	    
	    if ($f3 && /^Total Sequences\s+(\d+)/){

		$seq_tot = $1;

	    }
	    	
	    $f3 = 0 if /^>>END_MODULE/;
	    	
	}
    
        #translate low quality positions into intervals
        my $start = 0; 
        my $end = 0; 
        
        foreach my $pos (sort {$a <=> $b} @l_qual){
        
            if ($start == 0){
            
                $start = $pos;
                $end = $pos;
            
            } elsif ($pos == $end + 1){
            
                $end = $pos;
            
            } else {
		
		if ($start == $end){
				
		    $l_qual .= "$start ";
					
		} else {
				
                    $l_qual .= "$start-$end ";
		
		}
				
                $start = $pos;
		$end = $pos;

	    } #end of if 

	    if ($pos == $l_qual[-1]){
	
		if ($start == $end){

		    $l_qual .= "$start ";

		} else {

		    $l_qual .= "$start-$end ";

		}
			
	    } #end of if 
	    	
	} #end of foreach
	
	$sum{$sample}{$fastqc_report}{'0'} = "$seq_tot";
	$sum{$sample}{$fastqc_report}{'1'} = "$l_qual";
    
    }
    
}


#print extracted data into summary file

#open output file
open (OUT, ">$ms_report_dir/index.html");

print OUT "<HTML>";
print OUT "<TABLE><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Read name";
print OUT "<TH><CENTER>Number of reads";
print OUT "<TH><CENTER>Low quality positions";
print OUT "<TH><CENTER>Per base quality";
print OUT "<TH><CENTER>Per tile quality";
print OUT "<TH><CENTER>Per seq quality";
print OUT "<TH><CENTER>Per base content";
print OUT "<TH><CENTER>Per seq GC content";
print OUT "<TH><CENTER>Per base N content";
print OUT "<TH><CENTER>Seq length";
print OUT "<TH><CENTER>Seq duplication";
print OUT "<TH><CENTER>Overrepr seq";
print OUT "<TH><CENTER>Adapter content";
print OUT "<TH><CENTER>Kmer content\n";

my $f1 = 0;

#for each sample (alphanumerically sorted)
foreach my $sample (sort {$a cmp $b} keys %sum){

    #f1 is flagging for the first read group for each sample
    print OUT "<TR><TD>$sample";
    $f1 = 1;
    
    #for each read group
    foreach my $read (sort {$a cmp $b} keys %{$sum{$sample}}){

	#print error message
        if ($read eq "error"){

            if ($sum{$sample}{'error'} =~ /\w+/){

		print OUT "<TD COLSPAN=10>$sum{$sample}{'error'}";
		print OUT "\n";

	    }
	    
	    next;
	
	}
        
	#make a link to detailed fastqc report
        my $link = $summary_deployment;
        $link =~ s/www\/html\/(.*)$/$1/;
	$link = "http://"."$deployment_server/"."$link/"."$sample/$read/fastqc_report.html";
		
	#print read name
	my $read_title = $read;
	$read_title =~ s/^(\S+)_fastqc/$1/;
	
	if ($f1){
		
	    print OUT "<TD><A HREF = $link>$read_title</A>";
	    $f1 = 0;
		
	} else {
	    
	    print OUT "<TR><TD> <TD><A HREF = $link>$read_title</A>";
		
	}

	#print error message	
	if (defined($sum{$sample}{$read}{'error'})){
	    
	    print OUT "<TD COLSPAN=10>$sum{$sample}{$read}{'error'}";
	    print OUT "\n";
	    next;
		
	}
	
	#print FAIL/WARN/PASS for each quality metrics
	foreach my $cat (sort {$a <=> $b} keys %{$sum{$sample}{$read}}){
	
	    next unless $cat =~ /^\d+$/;
            
            if ($sum{$sample}{$read}{$cat} eq "FAIL"){
	        	
		print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";

            } elsif ($sum{$sample}{$read}{$cat} eq "WARN") {
            
            	print OUT "<TD><CENTER><IMG SRC=warning.png ALT=WARN>";
            	
            } elsif ($cat == 0||$cat == 1) {
	        
	        print OUT "<TD><CENTER>$sum{$sample}{$read}{$cat}</FONT>";
            
            } else {
	        
	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
	    
	    }

	}
		
	print OUT "\n";
    
    }

}

print OUT "</TABLE></HTML>";

#deploy HTML summary to the Web server
system("chmod 660 $ms_report_dir/index.html");
system("scp -r $ms_report_dir/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");


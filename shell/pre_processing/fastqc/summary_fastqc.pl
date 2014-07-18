#!/usr/bin/perl -w

$fastq_dir = "#pathReadsFastq";
$analysis_dir = "#pathAnalysisDir";
$report_dir = "#pathReportsDir";
$ms_report_dir = "#pathMSReportsDir";
$deployment_server = "#deploymentServer";
$summary_deployment = "#summaryDeployment";

#collect data from fastqc summaries in %sum
%sum = ();

#for each sample
opendir (PROJECT, "$fastq_dir")||print "Can't open $fastq_dir\n"; 
while (defined($sample = readdir(PROJECT))){
    next if $sample =~ /^\./;
    $sample_dir = "$fastq_dir/$sample";
    next unless (-d "$sample_dir");

    #check for FastQC run log file; integrity and proper pairing of fastq files
    if (-e "$analysis_dir/$sample/run/fastqc.$sample.log"){
        $sum{$sample}{'error'}=`grep ERROR $analysis_dir/$sample/run/fastqc.$sample.log`;
    }else{
	$sum{$sample}{'error'}="FastQC check was not done for this sample.";
    }

    next if $sum{$sample}{'error'} =~ /\w+/;

    #for each read retrieve FAIL/WARN/PASS for each quality metrics
    opendir (SAMPLE, "$sample_dir")||print "Can't open $sample_dir\n";
    while (defined($read = readdir(SAMPLE))){
        next unless $read =~ /^(\S+)\.f.*q.*/;
	$read = "$1"."_fastqc";
        if (-e "$report_dir/$sample/$read/summary.txt"){
	    $cat = 11;
	    open (SUM, "$report_dir/$sample/$read/summary.txt")||print "Can't open $report_dir/$sample/$read/summary.txt\n";
	    while (<SUM>){
                next if (/Basic Statistics/);
	        if (/^(\w\w\w\w)\t/){
		    $sum{$sample}{$read}{$cat} = $1;
	            $cat++;
		}
      	    }
	}else{
            $sum{$sample}{$read}{'error'}="FastQC check was not done for this read.";
	    next;
	}

        #read positions for low quality bases and total sequences produced for each read for each sample
	$f1 = 0; $f3 = 0;
        @l_qual = (); $l_qual = ""; 
        $seq_tot = 0;
	open (DATA, "$report_dir/$sample/$read/fastqc_data.txt")||print "Can't open $report_dir/$sample/$read/fastqc_data.txt\n";
	while (<DATA>){
	    $f1 = 1 if /^>>Per base sequence quality/;
	    if ($f1 && /^(\d+)\t\S+\t(\d+)\S+\t(\d+)\S+\t/){
                $pos = $1;
                $med = $2;
                $lq = $3;
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
        $start = 0; $end = 0; 
        foreach $pos (sort {$a <=> $b} @l_qual){
            if ($start == 0){
                $start = $pos;
                $end = $pos;
            }elsif ($pos == $end + 1){
                $end = $pos;
            }else{
		if ($start == $end){
                    $l_qual .= "$start ";
		}else{
                    $l_qual .= "$start-$end ";
		}
                $start = $pos;
		$end = $pos;
	    }
	    if ($pos == $l_qual[-1]){
		if ($start == $end){
                    $l_qual .= "$start ";
		}else{
                    $l_qual .= "$start-$end ";
		}
	    }
	}
	$sum{$sample}{$read}{'0'} = "$seq_tot";
	$sum{$sample}{$read}{'1'} = "$l_qual";
    }
}

#print extracted data into summary file
open (OUT, ">$ms_report_dir/index.html");
print OUT "<HTML>";
print OUT "<TABLE><TR>";
print OUT "<TH><CENTER>Sample";
print OUT "<TH><CENTER>Read name";
print OUT "<TH><CENTER>Number of reads";
print OUT "<TH><CENTER>Low quality positions";
print OUT "<TH><CENTER>Per base quality";
print OUT "<TH><CENTER>Per seq quality";
print OUT "<TH><CENTER>Per base content";
print OUT "<TH><CENTER>Per base GC content";
print OUT "<TH><CENTER>Per seq GC content";
print OUT "<TH><CENTER>Per base N content";
print OUT "<TH><CENTER>Seq length";
print OUT "<TH><CENTER>Seq duplication";
print OUT "<TH><CENTER>Overrepr seq";
print OUT "<TH><CENTER>Kmer content\n";
$f1 = 0;
foreach $sample (sort {$a cmp $b} keys %sum){
    print OUT "<TR><TD>$sample";
    $f1 = 1;
    foreach $read (sort {$a cmp $b} keys %{$sum{$sample}}){
        if ($read eq "error"){
            if ($sum{$sample}{'error'} =~ /\w+/){
	        print OUT "<TD COLSPAN=10>$sum{$sample}{'error'}";
                print OUT "\n";
	    }
	    next;
	}
        $link = $summary_deployment;
        $link =~ s/www\/html\/(.*)$/$1/;
	$link = "http://"."$deployment_server/"."$link/"."$sample/$read/fastqc_report.html";
	$read_title = $read;
	$read_title =~ s/^(\S+)_fastqc/$1/;
	if ($f1){
	    print OUT "<TD><A HREF = $link>$read_title</A>";
	    $f1 = 0;
	}else{
	    print OUT "<TR><TD> <TD><A HREF = $link>$read_title</A>";
	}
	if (defined($sum{$sample}{$read}{'error'})){
	    print OUT "<TD COLSPAN=10>$sum{$sample}{$read}{'error'}";
	    print OUT "\n";
	    next;
	}
	foreach $cat (sort {$a <=> $b} keys %{$sum{$sample}{$read}}){
	    next unless $cat =~ /^\d+$/;
            if ($sum{$sample}{$read}{$cat} eq "FAIL"){
	        print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
            }elsif ($sum{$sample}{$read}{$cat} eq "WARN"){
	        print OUT "<TD><CENTER><IMG SRC=warning.png ALT=WARN>";
            }elsif ($cat == 0||$cat == 1){
	        print OUT "<TD><CENTER>$sum{$sample}{$read}{$cat}</FONT>";
            }else{
	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
	    }
	}
	print OUT "\n";
    }
}
print OUT "</TABLE></HTML>";
system("chmod 660 $ms_report_dir/index.html");

system("scp -r $ms_report_dir/index.html $deployment_server:$summary_deployment/index.html > /dev/null 2>&1");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html > /dev/null 2>&1");



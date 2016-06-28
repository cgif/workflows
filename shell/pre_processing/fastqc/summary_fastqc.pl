#!/usr/bin/perl -w

use strict;
use warnings;

my $fastq_dir = "#pathReadsFastq";
my $runs_dir="#pathRunsDir";
my $report_dir = "#pathReportsDir";
my $ms_report_dir = "#pathMSReportsDir";
my $deployment_server = "#deploymentServer";
my $summary_deployment = "#summaryDeployment";

my @fields = split /\//, $runs_dir;
my $project_name = $fields[@fields-3];
my $run_date = $fields[@fields-1];

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
	my $run_dir = "$runs_dir/$sample";
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
				next if (/Basic Statistics/);
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
	    		} # end of if 
		} #end of foreach
	
		$sum{$sample}{$fastqc_report}{'0'} = "$seq_tot";
		$sum{$sample}{$fastqc_report}{'1'} = "$l_qual";
	} #end of while
} #end of main while



#print extracted data into summary file
open (OUT, ">$ms_report_dir/index.html");
print OUT"<html><head><title>FastQC Report</title><style type='text/css'> \@media screen {div.main {    display:block;    position:absolute;    overflow:auto;    height:auto;    width:auto;    top:4.5em;    bottom:2.3em;    left:0em;    right:0;    border-left: 1px solid #CCC;    padding:0 0 0 1em;    background-color: white;    z-index:1;  }    div.header {    background-color: #EEE;    border:0;    margin:0;    padding: 0.5em;    font-size: 200%;    font-weight: bold;    position:fixed;    width:100%;    top:0;    left:0;    z-index:2;  }  div.footer {    background-color: #EEE;    border:0;    margin:0;      padding:0.5em;    height: 1.3em;        overflow:hidden;    font-size: 100%;    font-weight: bold;    position:fixed;    bottom:0;    width:100%;    z-index:2;  }    img.indented {    margin-left: 3em;  } }  \@media print {  img {           max-width:100% !important;              page-break-inside: avoid;       }       h2, h3 {                page-break-after: avoid;        }       div.header {      background-color: #FFF;    }   }  body {      font-family: sans-serif;     color: #000;     background-color: #FFF;  border: 0;  margin: 0;  padding: 0;  }    div.header {  border:0;  margin:0;  padding: 0.5em;  font-size: 200%;  font-weight: bold;  width:100%;  }        #header_title {  display:inline-block;  float:left;  clear:left;  }  #header_filename {  display:inline-block;  float:right;  clear:right;  font-size: 50%;  margin-right:2em;  text-align: right;  }  div.header h3 {  font-size: 50%;  margin-bottom: 0;  } div.main {  background-color: white;  }        div.module {  padding-bottom:1.5em;  padding-top:1.5em;  }       div.footer {  background-color: #EEE;  border:0;  margin:0;  padding: 0.5em;  font-size: 100%;  font-weight: bold;  width:100%;  }  a {  color: #000080;  }  a:hover {  color: #800000;  }        h2 {  color: #800000;  padding-bottom: 0;  margin-bottom: 0;  clear:left;  }  table {   margin-left: 3em;  text-align: center;  }    th {   text-align: center;  background-color: #000080;  color: #FFF;  padding: 0.4em;  }          td {   font-family: monospace;   text-align: left;  background-color: #EEE;  color: #000;  padding: 0.4em;  }  img {  padding-top: 0;  margin-top: 0;  border-top: 0;  }    p {  padding-top: 0;  margin-top: 0;  }</style></head><body><div class='header'>       <div id='header_title'>         <img style='height:30px' src='igf.png' alt='IGF'/>FastQC Samples Report    </div>  <div id='header_filename'>".$run_date."<br>".$project_name."</div></div><div class='main'><div class='module'><h2 id='M0'></h2>";

print OUT "<table><thead><tr>";
print OUT "</tr></thead><tbody>";
print OUT "<th>MultiQC Report</th>";
print OUT "<th>FastQC Report</th>";
print OUT "</tr></thead><tbody>";
print OUT "<tr>";
print OUT "<td><center><a href=multiqc.html>MultiQC</a></FONT></td>";
print OUT "<td><center><a href=fastqc.html>FastQC</a></FONT></td>";
print OUT "</tr>";
print OUT "</tbody></table></div></div><div class='footer'></div></body></html>";

#open output file
open (OUT, ">$ms_report_dir/fastqc.html");

print OUT"<html><head><title>FastQC Report</title><style type='text/css'> \@media screen {div.main {    display:block;    position:absolute;    overflow:auto;    height:auto;    width:auto;    top:4.5em;    bottom:2.3em;    left:0em;    right:0;    border-left: 1px solid #CCC;    padding:0 0 0 1em;    background-color: white;    z-index:1;  }    div.header {    background-color: #EEE;    border:0;    margin:0;    padding: 0.5em;    font-size: 200%;    font-weight: bold;    position:fixed;    width:100%;    top:0;    left:0;    z-index:2;  }  div.footer {    background-color: #EEE;    border:0;    margin:0;      padding:0.5em;    height: 1.3em;        overflow:hidden;    font-size: 100%;    font-weight: bold;    position:fixed;    bottom:0;    width:100%;    z-index:2;  }    img.indented {    margin-left: 3em;  } }  \@media print {  img {           max-width:100% !important;              page-break-inside: avoid;       }       h2, h3 {                page-break-after: avoid;        }       div.header {      background-color: #FFF;    }   }  body {      font-family: sans-serif;     color: #000;     background-color: #FFF;  border: 0;  margin: 0;  padding: 0;  }    div.header {  border:0;  margin:0;  padding: 0.5em;  font-size: 200%;  font-weight: bold;  width:100%;  }        #header_title {  display:inline-block;  float:left;  clear:left;  }  #header_filename {  display:inline-block;  float:right;  clear:right;  font-size: 50%;  margin-right:2em;  text-align: right;  }  div.header h3 {  font-size: 50%;  margin-bottom: 0;  } div.main {  background-color: white;  }        div.module {  padding-bottom:1.5em;  padding-top:1.5em;  }       div.footer {  background-color: #EEE;  border:0;  margin:0;  padding: 0.5em;  font-size: 100%;  font-weight: bold;  width:100%;  }  a {  color: #000080;  }  a:hover {  color: #800000;  }        h2 {  color: #800000;  padding-bottom: 0;  margin-bottom: 0;  clear:left;  }  table {   margin-left: 3em;  text-align: center;  }    th {   text-align: center;  background-color: #000080;  color: #FFF;  padding: 0.4em;  }          td {   font-family: monospace;   text-align: left;  background-color: #EEE;  color: #000;  padding: 0.4em;  }  img {  padding-top: 0;  margin-top: 0;  border-top: 0;  }    p {  padding-top: 0;  margin-top: 0;  }</style></head><body><div class='header'>       <div id='header_title'>         <img style='height:30px' src='igf.png' alt='IGF'/>FastQC Samples Report    </div>  <div id='header_filename'>".$run_date."<br>".$project_name."</div></div><div class='main'><div class='module'><h2 id='M0'></h2>";
print OUT "<table><thead><tr>";
print OUT "<th>Sample</th>";
print OUT "<th>Read name</th>";
print OUT "<th>Number of reads</th>";
print OUT "<th>Low quality positions</th>";
print OUT "<th>Per base quality</th>";
print OUT "<th>Per tile quality</th>";
print OUT "<th>Per seq quality</th>";
print OUT "<th>Per base content</th>";
print OUT "<th>Per seq GC content</th>";
print OUT "<th>Per base N content</th>";
print OUT "<th>Seq length</th>";
print OUT "<th>Seq duplication</th>";
print OUT "<th>Overrepr seq</th>";
print OUT "<th>Adapter content</th>";
print OUT "<th>Kmer content</th>";
print OUT "</tr></thead><tbody>";
my $f1 = 0;

#for each sample (alphanumerically sorted)
foreach my $sample (sort {$a cmp $b} keys %sum){

    print OUT "<TR><TD>$sample</TD>";
    $f1 = 1;
    
    foreach my $read (sort {$a cmp $b} keys %{$sum{$sample}}){

        if ($read eq "error"){

            if ($sum{$sample}{'error'} =~ /\w+/){
				print OUT "<td COLSPAN=10>$sum{$sample}{'error'}</td>";
				print OUT "\n";
	    	}
	    	next;
	}
        
        my $link = $summary_deployment;
        $link =~ s/www\/html\/(.*)$/$1/;
	$link = "http://"."$deployment_server/"."$link/"."$sample/$read/fastqc_report.html";
		
	my $read_title = $read;
	$read_title =~ s/^(\S+)_fastqc/$1/;
	
	if ($f1){
	    	print OUT "<td><A HREF = $link>$read_title</A></td>";
	    	$f1 = 0;
	} else {
	    	print OUT "<TR><td></td><TD><A HREF = $link>$read_title</A></td>";
	}
	if (defined($sum{$sample}{$read}{'error'})){
	   	print OUT "<td COLSPAN=10>$sum{$sample}{$read}{'error'}</td>";
	    	print OUT "\n";
	    	next;
	}
	foreach my $cat (sort {$a <=> $b} keys %{$sum{$sample}{$read}}){
		next unless $cat =~ /^\d+$/;
		my $module=$cat-10;
		if ($sum{$sample}{$read}{$cat} eq "FAIL"){
	        	print OUT "<td><center><a href= $link"."#M".$module." ><IMG SRC=error.png ALT=FAIL></a></td>";
            	} elsif ($sum{$sample}{$read}{$cat} eq "WARN") {
            		print OUT "<td><center><a href= $link"."#M".$module."><IMG SRC=warning.png ALT=WARN></a></td>";
            	} elsif ($cat == 0||$cat == 1) {
	        	print OUT "<td><center>$sum{$sample}{$read}{$cat}</FONT></td>";
            	} else {
	        	print OUT "<td><center><a href= $link"."#M".$module."><IMG SRC=tick.png ALT=PASS></a></td>";
	    	}
	}
	print OUT "</TR>\n";
    }
}

print OUT "</tbody></table></div></div><div class='footer'></div></body></html>";

system("chmod -Rf 660 $ms_report_dir/*");
system("scp -r $ms_report_dir/* $deployment_server:$summary_deployment/");
system("ssh $deployment_server chmod -R 664 $summary_deployment/*");


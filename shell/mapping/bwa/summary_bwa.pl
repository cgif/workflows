#!/usr/bin/perl -w

$is_project_dir = "T";
$project_dir_analysis = "projectDirAnalysis";
$project_dir_results = "projectDirResults";
$today = "Today";
$summary_results = "summaryResults";
$deployment_server = "deploymentServer";
$summary_deployment = "summaryDeployment";

#collect data from bwa log files
%sum = (); %fastq = (); %mate = ();
if ($is_project_dir eq "T"){
$project_dir_analysis =~ s/^(\S+)\/\S+$/$1/;
$project_dir_results =~ s/^(\S+)\/\S+$/$1/;
opendir (PROJECT, "$project_dir_analysis"); 
while (defined($sample = readdir(PROJECT))){
    next if $sample =~ /^\./;
    $sample_dir = "$project_dir_analysis/$sample/run";
    next unless (-d "$sample_dir");

    if (`grep "does not contain any fastq files" $sample_dir/setup.log`){
        $fastq{$sample} = "FAIL";
    }else{
        $fastq{$sample} = "PASS";
    }

    if (`grep "No mate file found" $sample_dir/setup.log`){
        $mate{$sample} = "FAIL";
    }else{
        $mate{$sample} = "PASS";
    }

    opendir (SAMPLE, "$sample_dir"); 
    while (defined($log = readdir(SAMPLE))){
        next unless $log =~ /\.log/;
    
        if ($log =~ /^bwaAlignPe\.(\S+)\.f.*q\.(\w+)\.vs\..*\.log/){
	    $name = $1;
	    $split = $2;

            if (`grep "deleting temporary fastq files" $sample_dir/$log`){
                $sum{$sample}{$name}{$split}{'bam'} = "PASS";
            }else{
                $sum{$sample}{$name}{$split}{'bam'} = "FAIL";
            }

	}elsif ($log =~ /^samtoolsMerge\.((\S+)\.f.*q\.vs\..*)\.log/){
	    $full_name = $1;
	    $name = $2;

            if (`grep "deleting intermediate BAM files" $sample_dir/$log`){
                $sum{$sample}{$name}{'00'}{'merge'} = "PASS";
            }else{
                $sum{$sample}{$name}{'00'}{'merge'} = "FAIL";
            }

	    $stat = "$project_dir_results/$sample/$full_name".".sorted.bam.flagstat";
	    open (STAT,"$stat")||print "Can't open $stat\n";
	    while(<STAT>){
	        $sum{$sample}{$name}{'00'}{'total'} = "$1<BR>" if /^(\d+)\s+\+\s+\d+\s+in total/;
                $sum{$sample}{$name}{'00'}{'mapped'} = "$1($2)<BR>" if  /^(\d+)\s+\+\s+\d+\s+mapped \((\S+):.*/;
	    }
	}
    }
}
}else{
    $project_dir_analysis =~ /^\S+\/(\S+)\/\S+$/;
    $sample = $1;
    $sample_dir = "$project_dir_analysis/run";
    if (`grep "does not contain any fastq files" $sample_dir/setup.log`){
        $fastq{$sample} = "FAIL";
    }else{
        $fastq{$sample} = "PASS";
    }

    if (`grep "No mate file found" $sample_dir/setup.log`){
        $mate{$sample} = "FAIL";
    }else{
        $mate{$sample} = "PASS";
    }

    opendir (SAMPLE, "$sample_dir"); 
    while (defined($log = readdir(SAMPLE))){
        next unless $log =~ /\.log/;
    
        if ($log =~ /^bwaAlignPe\.(\S+)\.f.*q\.(\w+)\.vs\..*\.log/){
	    $name = $1;
	    $split = $2;

            if (`grep "deleting temporary fastq files" $sample_dir/$log`){
                $sum{$sample}{$name}{$split}{'bam'} = "PASS";
            }else{
                $sum{$sample}{$name}{$split}{'bam'} = "FAIL";
            }

	}elsif ($log =~ /^samtoolsMerge\.((\S+)\.f.*q\.vs\..*)\.log/){
	    $full_name = $1;
	    $name = $2;

            if (`grep "deleting intermediate BAM files" $sample_dir/$log`){
                $sum{$sample}{$name}{'00'}{'merge'} = "PASS";
            }else{
                $sum{$sample}{$name}{'00'}{'merge'} = "FAIL";
            }

	    $stat = "$project_dir_results/$full_name".".sorted.bam.flagstat";
	    open (STAT,"$stat");
	    while(<STAT>){
	        $sum{$sample}{$name}{'00'}{'total'} = "$1<BR>" if /^(\d+)\s+\+\s+\d+\s+in total/;
                $sum{$sample}{$name}{'00'}{'mapped'} = "$1($2)<BR>" if  /^(\d+)\s+\+\s+\d+\s+mapped \((\S+):.*/;
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
print OUT "<TH><CENTER>Fastq input";
print OUT "<TH><CENTER>Split";
print OUT "<TH><CENTER>BAM alignment";
print OUT "<TH><CENTER>Merging BAM";
print OUT "<TH><CENTER>Total reads";
print OUT "<TH><CENTER>Mapped reads\n";
$f1 = 0;
foreach $sample (sort {$a cmp $b} keys %sum){
    print OUT "<TR><TD>$sample";
    $f1 = 1;
    if ($fastq{$sample} eq "PASS"){
	print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
    }elsif ($fastq{$sample} eq "FAIL"){
	print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
    }

    if ($mate{$sample} eq "PASS"){
	print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
    }elsif ($mate{$sample} eq "FAIL"){
	print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
    }

    foreach $group (sort {$a cmp $b} keys %{$sum{$sample}}){
	if ($f1){
	    print OUT "<TD>$group";
	    $f1 = 0;
	}else{
	    print OUT "<TR><TD> <TD> <TD> <TD>$group";
	}

        foreach $split (sort {$a cmp $b} keys %{$sum{$sample}{$group}}){
	    if ($split eq '00'){
                print OUT "<TD>$split\n";
	    }else{
                print OUT "<TR><TD> <TD> <TD> <TD> <TD>$split\n";
	    }

            if ($sum{$sample}{$group}{$split}{'bam'} eq "PASS"){
	        print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
            }elsif ($sum{$sample}{$group}{$split}{'bam'} eq "FAIL"){
	        print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	    }
	 
            if ($split eq '00'){
                if ($sum{$sample}{$group}{$split}{'merge'} eq "PASS"){
	            print OUT "<TD><CENTER><IMG SRC=tick.png ALT=PASS>";
                }elsif ($sum{$sample}{$group}{$split}{'merge'} eq "FAIL"){
	            print OUT "<TD><CENTER><IMG SRC=error.png ALT=FAIL>";
	        }
                print OUT "<TD><CENTER>$sum{$sample}{$group}{$split}{'total'}" if defined($sum{$sample}{$group}{$split}{'total'});
                print OUT "<TD><CENTER>$sum{$sample}{$group}{$split}{'mapped'}" if defined($sum{$sample}{$group}{$split}{'mapped'});
            }

	    print OUT "\n";
	}
    }
}
print OUT "</TABLE></HTML>";

system("scp -r $summary_results/index.html $deployment_server:$summary_deployment/index.html");
system("ssh $deployment_server chmod -R 664 $summary_deployment/index.html");


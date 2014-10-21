#!/usr/bin/perl -w
use strict; 
use Getopt::Long;

my ($targets, $names, $HGMD_SNV, $HGMD_INDEL, $CLINVAR);
my $hgmdRelease;
my $clinvarDate;
my @geneList;
my ($out, $log, $cnvCalls);

$targets = "/groupvol/cgi/resources/annotations/Agilent_SureSelect_Human_All_Exon_V4+UTRs.targets.geneIDs.bed";
$names = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/GRCh37.p12.Ensemble_IDs.txt";

GetOptions ("hgmd=s" => \$hgmdRelease,
			"clinvar=s" => \$clinvarDate,
			"genes=s" => \@geneList,
			"out=s" => \$out,
			"log=s" => \$log,
			"input=s" => \$cnvCalls) or die("Error in command line arguments\n");

@geneList = split(/,/,join(',',@geneList)); #multiple gene list files are allowed as allow comma-separated lists of values as well as multiple occurrences of the option --genes

$hgmdRelease ||= '2014_3';
$clinvarDate ||= '2014_09_16';

$targets = "/groupvol/cgi/resources/annotations/Agilent_SureSelect_Human_All_Exon_V4+UTRs.targets.geneIDs.bed";
$names = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/GRCh37.p12.Ensemble_IDs.txt";
$HGMD_SNV = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/HGMD_Advanced_Substitutions.$hgmdRelease.csv";
$HGMD_INDEL = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/HGMD_Advanced_Micro_Lesions.$hgmdRelease.csv";
$CLINVAR = "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/clinvar_variant_summary_${clinvarDate}.txt";

open(OUT, ">$out") or die "Unable to open output file: $!";
open(LOG, ">$log") or die "Unable to open log file: $!";

my $header = `cat $cnvCalls|head -n 1 `;
chomp($header);
my @headerArray = split(/\t/, $header);
my $n = $#headerArray; 
print OUT "gene_name\tgene_ID\tpathogenic/associated_SNPs\tsample\ttype\tlocation\tBF\treads.ratio\t";
foreach my $annotation (@headerArray[14..$n]){
    print OUT "$annotation\t";
}
print OUT "\n";

my %count = ();
my %rev_count = ();
my %info_SNP = ();
my %info_CNV = ();
my %list = ();

foreach my $list (@geneList){

    open(LIST, "$list") or die "Unable to open gene list file $list: $!";
    while (my $gene = <LIST>){

        #check if $gene is an Ensembl id
	chomp($gene);
	unless ($gene =~ /ENSG/){
	    print LOG "Entry $gene is not an Ensembl id\n";
	    next;
	}

	#check if $gene is represented in exon targets 
	my $target = "";
	$target = `grep $gene $targets`;
	unless ($target ne ""){
	    print LOG "No exon targets for gene $gene\n";
	    next;
	}

	$list{$gene}++;
    }
}

foreach my $gene (keys %list){

    #extract gene name
    my $gene_name = ""; 
    $gene_name = `grep $gene $names|cut -f 3|sort|uniq`;
    chomp($gene_name);

    #collect data about pathogenic mutations
    my $grep = "";
	my @grep = ();
 
    $grep  = `grep $gene_name $HGMD_SNV`;
    @grep = split(/\n/, $grep);
    foreach my $line (@grep){
        my @info = split(/,/, $line);
	next if $info[1] eq "FP";
	$info_SNP{$gene} .= "$info[2] $info[9];";
    }

#    $grep = ""; @grep = (); 
    $grep  = `grep $gene_name $HGMD_INDEL`;
    @grep = split(/\n/, $grep);
    foreach my $line (@grep){
        my @info = split(/,/, $line);
	next if $info[1] eq "FP";
	$info_SNP{$gene} .= "$info[19] $info[3];";
    }

#    $grep = ""; @grep = (); 
    $grep  = `grep $gene_name $CLINVAR`;
    @grep = split(/\n/, $grep);
    foreach my $line (@grep){
        my @info = split(/\t/, $line);
		my @class = split(/;/,$info[5]);
		my $path = ""; 
		my $ass = "";
		foreach my $class (@class){
	    	$path++ if $class eq "pathogenic";
	    	$ass++ if $class eq "association";
		}    
		if ($path){
	    	$info_SNP{$gene} .= "$info[8] pathogenic $info[10];";
		}elsif ($ass){
	    	$info_SNP{$gene} .= "$info[8] association $info[10];";
		}
    }

#    $grep = ""; @grep = ();
    $grep = `grep $gene $cnvCalls`;
#	if ($grep ne ""){ 
    	@grep = split(/\n/, $grep);
    	foreach my $line (@grep){
			chomp $line;
			my @info = split(/\t/, $line);
        	my $sample = $info[0];
        	$info_CNV{$gene}{$gene_name}{$sample} .= "$info[2]\t$info[1]\t$info[8]\t$info[11]\t";
			foreach my $annotation (@info[14..$n]){
	     		$info_CNV{$gene}{$gene_name}{$sample} .= "$annotation\t";
			}
			$info_CNV{$gene}{$gene_name}{$sample} .= "|";
			$count{$gene}++;
    	}
#	}
}

foreach my $gene (keys %count){
    $rev_count{$count{$gene}} .= "$gene|";
}

foreach my $count (sort { $a <=> $b } keys %rev_count){
    my @genes = split(/\|/, $rev_count{$count});
    foreach my $gene (@genes){
	next if $gene eq "";
        foreach my $gene_name (keys %{$info_CNV{$gene}}){ 
	    	if (defined($info_SNP{$gene})){
        		print OUT "$gene_name\t$gene\t$info_SNP{$gene}\t";
	    	}else{
        		print OUT "$gene_name\t$gene\tNA\t";
	    	}
	    	my $line = 0;
        	foreach my $sample (keys %{$info_CNV{$gene}{$gene_name}}){ 
				my @CNVs = split(/\|/, $info_CNV{$gene}{$gene_name}{$sample});
				foreach my $CNV (@CNVs){
		    		next if $CNV eq "";
		    		unless ($line){
						print OUT "$sample\t$CNV\n";
		    		}else{
						print OUT " \t \t \t$sample\t$CNV\n";
		    		}
		    		$line++;
				}
	    	}
		}
	}
}

close (OUT) or die "Unable to close output file: $!";
close (LOG) or die "Unable to close log file: $!";


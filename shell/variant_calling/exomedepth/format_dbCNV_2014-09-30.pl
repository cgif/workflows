#!/usr/bin/perl -w

#open (DATA, "/groupvol/cgi/resources/annotations/GRCh37_hg19_DGV_variants.2013-07-23.txt")||print "Can't open GRCh37_hg19_DGV_variants.2013-07-23.txt\n";
#open (BED, ">/groupvol/cgi/resources/annotations/GRCh37_hg19_DGV_variants.2013-07-23.multisampleCNV.bed");
#while(<DATA>){
#    chomp();
#    @data = split(/\t/,$_);
#    next unless "$data[4]" eq "CNV";
#    @length = split(/,/,$data[19]);
#    $length = @length;
#    if ($data[12] eq 'M' || $data[15] > 1 || $data[16] > 1 || @length > 1){
#        print BED "$data[1]\t$data[2]\t$data[3]\t$data[0]\n";
#    }
#}

######################

#open (DATA_CUR, "/groupvol/cgi/resources/annotations/ISCA_Curated_Pathogenic_CNVs_hg19_2011-08-25.bed")||print "Can't open ISCA_Curated_Pathogenic_CNVs_hg19_2011-08-25.bed\n";
#open (DATA_PAT, "/groupvol/cgi/resources/annotations/ISCA_Pathogenic_CNVs_hg19_2012-06-14.bed")||print "ISCA_Pathogenic_CNVs_hg19_2012-06-14.bed\n";

#%data = ();
#while(<DATA_CUR>){
#    chomp();
#    next unless /^chr/;
#    @data = split(/\t/,$_);
#    $chrom = $data[0];
#    $chrom =~ s/chr(.*)/$1/;
#    $data{$chrom}{$data[1]}{$data[3]} = "$chrom\t$data[1]\t$data[2]\t$data[3]\n";
#}

#while(<DATA_PAT>){
#    chomp();
#    next unless /^chr/;
#    @data = split(/\t/,$_);
#    $chrom = $data[0];
#    $chrom =~ s/chr(.*)/$1/;
#    $data{$chrom}{$data[1]}{$data[3]} = "$chrom\t$data[1]\t$data[2]\t$data[3]\n";
#}

#open (BED, ">/groupvol/cgi/resources/annotations/ISCA_Pathogenic+Curated_CNVs_hg19_2012-06-14.bed");
#foreach $chrom (qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)){
#    foreach $start (sort {$a <=> $b} keys %{$data{$chrom}}){
#        foreach $name (keys %{$data{$chrom}{$start}}){
#	    print BED "$data{$chrom}{$start}{$name}";
#        }
#    }
#}

##################

open (DATA, "/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/clinvar_variant_summary_2014_09_16.txt")||print "Can't open clinvar_variant_summary_2014_09_16.txt\n";
open (BED, ">/groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/clinvar_variant_summary_2014_09_16.pathogenic.50bp.bed");

while(<DATA>){
    chomp();
    @data = split(/\t/,$_);

    @cs = split(/;/,$data[5]);
    $pathogenic = 0;
    foreach $cs (@cs){
	$pathogenic++ if "$cs" eq "Pathogenic";
    }
    next unless $pathogenic;

    next unless $data[15] > 1;
    $length = $data[15] - $data[14];
    next unless $length > 50;    

    print BED "$data[13]\t$data[14]\t$data[15]\t$data[8]\n";
}

###################

#open(DATA,"/home/asoskins/projects/blakemore_obesitywe/ISCA.obesity.csv");
#open(BED,">/home/asoskins/projects/blakemore_obesitywe/ISCA.obesity.bed");
#while(<DATA>){
#    chomp();
#    next unless /^nssv/;
#    @data = split(/\t/,$_);
#    $loci = "$data[1]";
#    $inter = "pathog" if "$data[7]" eq "Pathogenic";
#    $inter = "benign" if "$data[7]" eq "Benign";
#    $inter = "uncert" if $data[7] =~ /Uncertain/;
#    $name = "$data[0]_$inter:$data[3]";
#    $loci =~ /chr(.*)\:(\d+)-(\d+)/;
#    print BED "$1\t$2\t$3\t$name\n";
#}

###################

#awk '{if (/#/) {next} if ($5 > 0) {print $4 "\t" $2 "\t" $3 "\t" $1 ":gain"} if ($5 < 0) {print $4 "\t" $2 "\t" $3 "\t" $1 ":loss"}}' /groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/decipher-hg19_14-09-18.txt > /groupvol/cgi/resources/annotations/eukaryote/human/GRCh37/CNV_annotations/decipher-hg19_14-09-18.bed

#grep obesity decipher-hg19_13-11-26.txt|awk '{if (/#/) {next} if ($5 > 0) {print $4 "\t" $2 "\t" $3 "\t" $1 ":gain"} if ($5 < 0) {print $4 "\t" $2 "\t" $3 "\t" $1 ":loss"}}' > /home/asoskins/projects/blakemore_obesitywe/obesity_targets/decipher-hg19_13-11-26.obesity.bed

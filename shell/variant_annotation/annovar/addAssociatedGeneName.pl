#!/usr/bin/perl

#reads a tab delimited file in the format:
#Ensembl Gene ID\tAssociated Gene Name
#and adds gene names as an additional column
#with the gene name to Ensembl based variant
#annotations

my $annovar_file=$ARGV[0];
my $gene_names_file=$ARGV[1];

#read associated gene names
my %gene_names;

my $gene_names_input;
if($gene_names_file =~ /\.gz/){
    open($gene_names_input, "gzip -cd $gene_names_file |")
}  else {
    open($gene_names_input,  "<$gene_names_file");
}

while(<$gene_names_input>){

    chomp;
    my @cols=split(/\t/);
    $gene_names{$cols[0]} = $cols[1];
    
}

#foreach my $key (keys %gene_names){
#   print "$key\t$gene_names{$key}\n";
#}

my $annovar_input;
if($annovar_file =~ /\.gz/){
    open($annovar_input, "gzip -cd $annovar_file |") || die "Unable to open Annovar file: $!";
}  else {
    open($annovar_input,  "<$annovar_file") || die "Unable to open Annovar file: $!";
}

while(<$annovar_input>){

    chomp;
    my @cols=split(/\t/);
    my $col_count=@cols;
    my @ensembl_ids=split(/[,;]/,$cols[6]);
    my $id_count=@ensembl_ids;
    my $names="";

    for (my $i=0; $i < $id_count; $i++){
	my $id = $ensembl_ids[$i];
	my $name = $gene_names{$id};
	if($name eq ""){
	    $name="NULL";
	}
	$names=$names.$name;
	if($i != $id_count-1){
	    $names=$names.",";
	}
    }
    
    #print first 6 columns
    for(my $i=0;  $i < 7;  $i++){
	print "$cols[$i]\t";
    }
    
    #insert name columns
    print "$names\n";

    #print remaining columns
     for (my $i=7; $i < $col_count; $i++){
	 print "$cols[$i]";
	if($i != $col_count-1){
	    print "\t";
	} else {
	    print "\n";
	}
    }


}


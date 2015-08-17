#!/usr/bin/perl -w

# Script for translation of gff file into refFlat (genePred that associates the gene name with the gene prediction information)
# In alternative splicing situations each transcript has a row in this table with the same gene name in the first column
# UCSC utility 'ldHgGene' generate genePred file with numbers insted of gene names in the first column
# Proper gene names are extracted from input gff file

$tmp_dir = $ARGV[0];
$gff = "$tmp_dir/tmp.gff";
$gp = "$tmp_dir/tmp.genePred";
$rf = "$tmp_dir/tmp.refFlat";

open(GP, "$gp");
open(RF, ">$rf");
while(<GP>){
    chomp();
    /^\S+(\s+(\S+)\s+.*)/;
    $info = $1; $tr_id = $2;
    $gene_name = "";
    $gene_name = `grep $tr_id $gff|head -n 1|cut -f 9`;
    $gene_name =~ s/.*gene_name \"\s*(\S+)\s*\".*/$1/;
    chomp($gene_name);
    print "$info\n" unless $gene_name;
    print RF "$gene_name$info\n";
}

#!/usr/bin/perl

use strict;
use warnings;

use lib '/project/tgu/libs/perl/lib/perl5';
#use Spreadsheet::WriteExcel;
#use Spreadsheet::WriteExcelXML; 

#use XLSX writer. The main advantage of the XLSX format over the XLS format 
#is that it allows a larger number of rows and columns in a worksheet. The 
#XLSX file format also produces much smaller files than the XLS file format. 
use Excel::Writer::XLSX;

my $tsv = $ARGV[0];
my $xls = $ARGV[1];

#no longer required if we use the XLSX format
#my $max_rows_xls=65536;
#
#open(TSV, "<$tsv");
#my $row_count = 0;
#while(<TSV>){
#	$row_count++;
#}
#close(TSV);

#if($row_count > $max_rows_xls){ 
#	$workbook = Spreadsheet::WriteExcelXML->new($xls);
#} else {
#	$workbook = Spreadsheet::WriteExcel->new($xls);
#}


# Create a new Excel workbook
my $workbook;

$workbook = Excel::Writer::XLSX->new($xls);

# Add a worksheet
my $worksheet = $workbook->add_worksheet();

#  Add and define a format
my $format_header = $workbook->add_format(); # Add a format
$format_header->set_bold();

# write tsv to xls
open(TSV, "<$tsv");

my $row = 0;
while(<TSV>){

	chomp;
	my @cols = split(/\t/);

	for(my $col = 0; $col < @cols; $col++){
		
		my $value = $cols[$col];
		
		#if header
		if($row == 0){		
			$worksheet->write($row, $col, $value, $format_header);
		} 
		else {
			$worksheet->write($row, $col, $value);			
		}
								
	}

	$row++;

}

close(TSV);


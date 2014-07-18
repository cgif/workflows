<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  
	<title>sample_summary.php</title>
  
	<style>
		
		h1 {
			font-family:"Arial", Arial, sans-serif;
		}
	
		table.custom {	
			font-family:"Arial", Arial, sans-serif;
			border-width: 1px;
			border-style: solid;
			border-color: #3399FF;
			border-collapse: collapse;
			text-align:center;
			
		}
		table.custom th {
			background-color: #3399FF;
			color: white;
			padding: 4px;
		}
		table.custom tr {
			color: black;
			padding: 4px;
		}
		tr.customOdd {
			background-color: #E8E8E8;		
		}
		
	</style>
	
</head>

<body>

<h1>#header</h1>
<div>
<?php

$data = file("#tsvFile");

print "<table class=\"custom\">";

#construct table
$m=0;
foreach ($data as $line){

	$m++;
	$lineArray = explode("\t", $line);

	#if header
	if($m == 1){
   
		print "<tr>";
		foreach($lineArray as $value){
			print "<th>$value</th>";
		}
		print "</tr>\n";

	}
	#if values
	else {

		#alternate row color
		if($m % 2){
			print "<tr class=\"customOdd\">";
		}
		else {
			print "<tr>";

		}

		foreach($lineArray as $value){
						
			print "<td>$value</td>";
			
		}
		print "</tr>\n";
	
	}
	
} // end foreach

//print the bottom of the table
print "</table> \n";

?>
</div>
    
</body>
</html>


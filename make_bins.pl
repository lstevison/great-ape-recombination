#! /usr/bin/perl

#May 10, 2013
#Modified September 22, 2014
#Program reads in BED file and thins blocks greater than user-defined scale into smaller blocks to perform comparative map analysis
#this input should be used to compute mean rate across reduced map files in user-defined scale windows

$input = $ARGV[0];
$scale=$ARGV[1];
$label=$ARGV[2];
$threshold=$scale*.9; #define threshold based on 90% coverage
unless ($#ARGV==2) {
    print STDERR "Please provide on command line an input BED file, a user defined scale and a label for output\n\n";
    die;
} #end unless

#$output = 'Gorilla_map_splitby' . $label . '_chr2b.BED';
$output = $input . ".binned_by_" . $label . ".txt";

open(INPUT, $input);
open(OUTPUT, ">$output");

print STDERR "Reading input...";

while(<INPUT>) {
    chomp;
    @input_line = split(/\t/, $_);
    $block_size = $input_line[2]-$input_line[1]+1;
    $loop_size = int($block_size/$scale)+1;
    $adjusted_start=1+int($input_line[1]/($scale+1))*$scale;
    $adjusted_end=($scale+1)+(int($input_line[2]/($scale+1))*$scale);

#    print STDERR "Adjusted start: $adjusted_start, $adjusted_end, $loop_size\n";

    for ($i=0; $i<=$loop_size; $i++) {
	$count=$i+1;
	$scaling_term=$scale*($count);
	
	if ($i==0) {
#	    $start=1+int($input_line[1]/($scale+1))*$scale;
	    $start=$adjusted_start;
	    $end=$adjusted_start+$scaling_term-1;
	    $last_end=$end;
	} elsif ($i==$loop_size) {
	    $start=$last_end+1;
	    $end=$adjusted_end;
	} else {
	    $start=$last_end+1;
	    $end=$adjusted_start+$scaling_term-1;
	    $last_end=$end;
	} #end else
	
	$region_size = $end-$start+1;
	if ($region_size >=$threshold) {
	    print OUTPUT "$input_line[0]\t$start\t$end\t$region_size\t$input_line[3]\t$input_line[4]\t$input_line[5]\n";
	} #end if
	
    } #end for

} #end while

print STDERR "done.\n";


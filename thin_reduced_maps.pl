#! /usr/bin/perl

#May 10, 2013
#Modified Apr 2, 2014
#Program reads in 'intersected_above1Mb_plusOri_fixed.BED' and thins blocks greater than user-defined scale into smaller blocks to perform comparative map analysis
#this input should be used to compute mean rate across reduced map files in user-defined scale windows

$input = 'intersected_above1Mb_plusOri_fixed.BED';
$scale=$ARGV[0];
$label=$ARGV[1];
$threshold=$scale*.9; #define threshold based on 90% coverage
unless ($#ARGV==1) {
    print STDERR "Please provide on command line a user defined scale and a label for output\n\n";
    die;
} #end unless

$output = 'intersected_above1Mb_fixed_plusOri_splitby' . $label . '.BED';

open(INPUT, $input);
open(OUTPUT, ">$output");

print STDERR "Reading input...";

while(<INPUT>) {
    chomp;
    @input_line = split(/\t/, $_);
    $block_size = $input_line[2]-$input_line[1]+1;
    $loop_size = int($block_size/$scale);

    for ($i=0; $i<=$loop_size; $i++) {
	$count=$i+1;
	$scaling_term=$scale*($count);
	
	if ($i==0) {
	    $start=$input_line[1];
	    $end=$input_line[1]+$scaling_term-1;
	    $last_end=$end;
	} elsif ($i==$loop_size) {
	    $start=$last_end+1;
	    $end=$input_line[2];
	} else {
	    $start=$last_end+1;
	    $end=$input_line[1]+$scaling_term-1;
	    $last_end=$end;
	} #end else
	
	$region_size = $end-$start+1;
	if ($region_size >=$threshold) {
	    print OUTPUT "$input_line[0]\t$start\t$end\t$region_size\t$input_line[3]\t$input_line[4]\t$input_line[5]\n";
	} #end if
	
    } #end for

} #end while

print STDERR "done.\n";


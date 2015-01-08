#! /usr/bin/perl

#program to match fimo output with input of matched ohtspots and coldspots

#Created by Laurie Stevison
#Created on: January 7, 2014
#Last modified: September 10, 2014

$matched=$ARGV[0]; #region coordinates match with rate and GC estimates
$fimo1=$ARGV[1]; #fimo output
$output=$ARGV[2]; #output filename
$output2=$output;
$output2=~s/.txt/.BED/;

unless ($#ARGV==2) {
    print STDERR "Please provide on command line: (1) extended BED file with rate and GC info, (2) fimo output, and (3) output filename.\n\n";
    die;
} #end unless

print STDERR "Reading in fimo outputs...";

open(FIMO_HOT, $fimo1);
$header=(<FIMO_HOT>);

@hotspot_name=();
%hotspot_pbs=();

while (<FIMO_HOT>) {
    chomp;
    @input_array2=split(/\s+/,$_);
    push(@hotspot_name, $input_array2[1]);
#    $mid=($input_array2[2]+$input_array2[3])/2;

    push @{$hotspot_pbs{$input_array2[1]}}, "$input_array2[2]:$input_array2[3]";
} #end while

open(MATCHES, $matched);
open(OUTPUT, ">$output");
open(OUTPUT2, ">$output2");

print STDERR "done.\nNow reading in regions matched with recombination and GC info and compiling fimo data into single output...";
print OUTPUT "chr\tstart\tend\tpeakRate\tavgRate\t%GC\t%N\tNum_hits\n";
print OUTPUT2 "chr\tstart\tend\n";

$header=(<MATCHES>);

while (<MATCHES>) {
    chomp;
    @input_array1=split(/\s+/,$_);
    $hotstart = $input_array1[1];
    $hotend=$input_array1[2];
#    print STDERR "$hotstart, $hotend, $hotname\n";
    $hotname = "$input_array1[0]\.$hotstart";
    @hot_hits = @{$hotspot_pbs{$hotname}}; #extracts fimo hits for hotspot

    if (defined(@hot_hits)) {
	$num_hothits = $#hot_hits+1;
	for ($j=0; $j<=$#hot_hits;$j++) {
	    @pbsite=split(":",$hot_hits[$j]);
	    $pbstart=$hotstart+$pbsite[0];
	    $pbend=$hotstart+$pbsite[1];
	    print OUTPUT2 "$input_array1[0]\t$pbstart\t$pbend\n";
	} #end for
    } else {
	$num_hothits = 0;
    } #end if

    print OUTPUT "$input_array1[0]\t$hotstart\t$hotend\t$input_array1[3]\t$input_array1[4]\t$input_array1[7]\t$input_array1[8]\t$num_hothits\n";

} #end while

print STDERR "done.\n";

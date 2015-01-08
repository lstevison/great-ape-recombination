#! /usr/bin/perl

#program to match fimo output with input of matched ohtspots and coldspots

#Created by Laurie Stevison
#Created on: January 7, 2014
#Last modified: August 28, 2014

$matched=$ARGV[0]; #matched coldspots and hotspots file from panmap
$fimo1=$ARGV[1]; #fimo output for hotspots
$fimo2=$ARGV[2]; #fimo output for coldspots
$output=$ARGV[3]; #output filename

unless ($#ARGV==3) {
    print STDERR "Please provide on command line: (1) matched hotspots and coldspots, (2) fimo output for hotspots, (3) fimo output for coldspots, and (4) output filename.\n\n";
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
    $mid=($input_array2[2]+$input_array2[3])/2;
    push @{$hotspot_pbs{$input_array2[1]}}, $mid;
} #end while

open(FIMO_COLD, $fimo2);
$header=(<FIMO_COLD>);

@coldspot_name=();
%coldspot_pbs=();

while (<FIMO_COLD>) {
    chomp;
    @input_array3=split(/\s+/,$_);
    push(@coldspot_name, $input_array3[1]);
    $mid=($input_array3[2]+$input_array3[3])/2;
    push @{$coldspot_pbs{$input_array3[1]}}, $mid;
} #end while

open(MATCHES, $matched);
open(OUTPUT, ">$output");

print STDERR "done.\nNow reading in matched hotspots and coldspots and compiling data into single output...";
print OUTPUT "Hotspot_chr\tHotspot_start\tHotspot_end\tHotspot_GC\tHotspot_N\tHotspot_avgRate\tHotspot_peakRate\tNum_hot_hits\tColdspot_chr\tColdspot_start\tColdspot_end\tColdspot_GC\tColdspot_N\tColdspot_avgRate\tColdspot_peakRate\tNum_cold_hits\n";

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
    } else {
	$num_hothits = 0;
    } #end if

    $coldstart=$input_array1[8];
    $coldend=$input_array1[9];
    $coldname = "$input_array1[0]\.$coldstart";
    @cold_hits = @{$coldspot_pbs{$coldname}}; #extracts fimo hits for coldspot

    if (defined(@cold_hits)) {
	$num_coldhits = $#cold_hits+1;
    } else {
	$num_coldhits = 0;
    } #end if

    print OUTPUT "$input_array1[0]\t$hotstart\t$hotend\t$input_array1[3]\t$input_array1[4]\t$input_array1[5]\t$input_array1[6]\t$num_hothits\t$input_array1[0]\t$coldstart\t$coldend\t$input_array1[10]\t$input_array1[11]\t$input_array1[12]\t$input_array1[13]\t$num_coldhits\n";

} #end while

print STDERR "done.\n";

#! /usr/bin/perl

#Program reads in get-INFO output from vcftools and kept.sites output from
#Gorilla final recode
#and re-prints Gorilla coords matched with hg18 coordinates

#Modified November 7, 2013

$input_INFO = $ARGV[0];
$input_Gorilla = $ARGV[1];
$output = $ARGV[2];

unless ($#ARGV==2) {
    print STDERR "Please provide on command line input INFO file from vcftools, Gorilla kept.sites output, and output filename\n\n";
    die;
} #end unless

open(INFO, $input_INFO);
open(GOR, $input_Gorilla);
open(OUTPUT, ">$output");

@hg18_start = ();
@hg18_chr = ();
@gorgor_start = ();

unless (-e $input_INFO) {
    print STDERR "File $input_INFO doesn't Exist!\n";
    die;
} #end unless

unless (-e $input_Gorilla) {
    print STDERR "File $input_Gorilla doesn't Exist!\n";
    die;
} #end unless

$header1 = (<INFO>);  #discard header of INFO
$header2 = (<GOR>);   #discard header of kept.sites

print STDERR "Now reading in INFO input file...";

while (<INFO>) {
    chomp;
    @line_item = split(/\t/, $_);
    $gorilla = $line_item[1]/1;

    if ($line_item[4]=~/"failed"/) {
	next;
    } else {
	@dummy = split(":", $line_item[4]);
        push(@hg18_chr, $dummy[0]);
        push(@hg18_start, $dummy[1]);
	push(@gorgor_start, $gorilla); 
    } #end else
} #end while

$i=0;
print STDERR "no of lines: $#hg18_chr, $#hg18_start, $#gorgor_start  ";
print STDERR "done.\nNow reading in kept sites output and writing final output file...";

while(<GOR>) {
    chomp;
    @gor_line = split(/\s+/, $_);
    $chr = $gor_line[0];
    if ($i==0) {
	print STDERR "$chr";
    } #end if

    for ($j=$i; $j<=$#hg18_start; $j++) {
	if ($gor_line[1]==$gorgor_start[$j]) {
	    print OUTPUT "$chr\t$gor_line[1]\t$hg18_chr[$j]\t$hg18_start[$j]\n";
	    $i = $j+1;
	    last;
	} #end if
    } #end for
    
} #end while

print "done.\n\n";


#! /usr/bin/perl

$burn_in = $ARGV[0];
$input1 = $ARGV[1];
$input2 = $ARGV[2];
$output = $ARGV[3];


#kill program if not enough command arguments
unless ($#ARGV==3) {
    print STDERR "Must specify output name on command line\n";
    die;
} #end unless

# open inputs
open(INPUT1, $input1);
open(INPUT2, $input2);

#dump header line of both files
$header1 = (<INPUT1>);
$header2 = (<INPUT2>);

#read in first input
while (<INPUT1>) {
    chomp;
    $line_item1 = $_;
    @line1 = split(/\s+/, $line_item1);
    push(@input1_array, $line1[0]);
} #end while

#read in second input
while (<INPUT2>) {
    chomp;
    @line2 = split(/\s+/,$_);
    push(@input2_array, $line2[0]);
} #end while

#print output
open(OUTPUT,">$output");

#print header of output
#print OUTPUT "total map length\t number of blocks\n";

for ($i=$burn_in; $i<=$#input1_array; $i++) {
    print OUTPUT "$input1_array[$i]\t$input2_array[$i]\n";
} #end for

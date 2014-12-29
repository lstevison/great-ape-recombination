#! /usr/bin/perl

#February 16, 2012
#Modified August 15, 2012

#program reads in BED file with Gorilla coordinates, hg18 coordinates
#and recently lifted-over reciprocal liftover results for hg18Togorilla

#program will calculate distance between new gorilla coordinates

$input = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
	print STDERR "Error: please provide input and output filenames on command line\n";
	die;
} #end unless

open (INPUT, $input);
open (OUTPUT, ">$output");

#print header
print OUTPUT "chromosome\tRL_gor_start\tRL_gor_end\tRL_dist\thg18_chr\thg18_start\thg18_end\thg18_dist\tORG_gor_chr\tORG_gor_start\tORG_gor_end\tORG_gor_dist\n";

while (<INPUT>) {
	chomp;
	@line_item = split(/\t/, $_);
	$new_gordist = $line_item[2] - $line_item[1];
	print OUTPUT "$line_item[0]\t$line_item[1]\t$line_item[2]\t$new_gordist\t$line_item[3]\t$line_item[4]\t$line_item[5]\t$line_item[6]\t$line_item[7]\t$line_item[8]\t$line_item[9]\t$line_item[10]\n";
} #end while

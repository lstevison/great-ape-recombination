#! /usr/bin/perl

#February 16, 2012
#Modified August 15, 2012

#program reads in BED file with original Gorilla coordinates and recently 
#lifted-over coordinates to hg18

#program will copy hg18 coordinates and put the 2nd copy at beginning
#for subsequent reciprocal liftover (retaining original coordinates)
#program will add two new columns with the distance between coordinates 
#in gorilla vcf file and hg18 coordinates

$input = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
	print STDERR "Error: please provide input and output filenames on command line\n";
	die;
} #end unless

print STDERR "Running liftover cleanup program...\n";

open (INPUT, $input);
open (OUTPUT, ">$output");

while (<INPUT>) {
	chomp;
	@line_item = split(/\t/, $_);
	$gor_dist = $line_item[2] - $line_item[1];
	$human_dist = $line_item[5] - $line_item[4];
	print OUTPUT "$line_item[0]\t$line_item[1]\t$line_item[2]\t$line_item[0]\t$line_item[1]\t$line_item[2]\t$human_dist\t$line_item[3]\t$line_item[4]\t$line_item[5]\t$gor_dist\n";
} #end while

print STDERR "\tdone.\n";

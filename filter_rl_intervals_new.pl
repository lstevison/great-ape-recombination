#! /usr/bin/perl

#February 16, 2012
#Modified August 15, 2012

#program reads in BED file with Gorilla coordinates, hg18 coordinates
#and recently lifted-over reciprocal liftover results for hg18Togorilla

#program will calculate distance between new gorilla coordinates

$input = $ARGV[0];
$original_input = $ARGV[1];
$output = $ARGV[2];

unless ($#ARGV==2) {
	print STDERR "Error: please provide input, original input, and output filenames on command line\n";
	die;
} #end unless

print STDERR "Running program to filter reciprocal liftover sites...";

open (INPUT, $input);
open(BED, $original_input);
open (OUTPUT, ">$output");
open(CLEAN, ">$output.include.txt");

$header = (<INPUT>);    #discard header line
@sites_to_keep = ();
print CLEAN "Chr\tposition\n";    #prints header to output

while (<INPUT>) {
	chomp;
	@line_item = split(/\t/, $_);
	$dont_filter_out=0;
	$chrom_match = 0;
	$start_match = 0;
	$end_match = 0;
	$length_match = 0;

	if ($line_item[0] eq $line_item[8]) {    #chr field matches between original gorilla and new gorilla
	    #$dont_filter_out = 1;
	    $chrom_match = 1;
	} else { 
	    next;
	} # end if

	#is original start within 100bp of new start
	#$lower_limit_start = $line_item[1] - 100;
	#$upper_limit_start = $line_item[1] + 100;

	if ($chrom_match==1 && $line_item[1]==$line_item[9]) {  #matches exactly original start
	    #$dont_filter_out=1;
	    $start_match = 1;
	} else {
	    next;
	} #end else

	#is original end within 100bp of new end
	#$lower_limit_end = $line_item[2] - 100;
	#$upper_limit_end = $line_item[2] + 100;

	if ($start_match==1 && $line_item[2]==$line_item[10]) {  #matches exactly original start
	    #$dont_filter_out=1;
	    $end_match = 1;
	} else {
	    next;
	} #end else

	#is gorgor length within 20bp of human length
	$lower_limit_length1 = $line_item[3] - 20;
	$upper_limit_length1 = $line_item[3] + 20;

	if ($end_match==1 && $line_item[7] < $upper_limit_length1 && $line_item[7] > $lower_limit_length1) {
	    #$dont_filter_out=1;
	    $length_match=1;
	} else {
	    next;
	} #end else

	#is gorgor original length within 20bp of new length
	$lower_limit_length2 = $line_item[3] - 20;
	$upper_limit_length2 = $line_item[3] + 20;

	if ($length_match==1 && $line_item[11] < $upper_limit_length2 && $line_item[11] > $lower_limit_length2) {
	    $dont_filter_out=1;
	} else {
	    next;
	} #end else
	
	if ($dont_filter_out==1) {
#	    print OUTPUT "$line_item[8]\t$line_item[9]\n";
	    $chromosome = $line_item[0];
	    $interval_start = "$chromosome:$line_item[9]";
#	    $interval_end = "$chromosome:$line_item[10]";
	    push (@sites_to_keep, $interval_start);
#	    push (@sites_to_keep, $interval_end);
	} #end if

} #end while

@sites_to_keep2 = ();

for ($z=0; $z<=$#sites_to_keep; $z++) {
    @site = split(":",$sites_to_keep[$z]);

    if ($site[0] eq $chromosome) {
	push(@sites_to_keep2, $site[1]);
    } #end if
} #end for

@sites_to_keep2 = sort { $a <=> $b } (@sites_to_keep2);
%hash = map { $_ => 1 } @sites_to_keep2;
@unique_sites = keys %hash;
@unique_sites = sort { $a <=> $b } (@unique_sites);

for ($k=0; $k<=$#unique_sites; $k++) {
    print CLEAN "$chromosome\t$unique_sites[$k]\n";
} #end for

$i=0;
    
while (<BED>) {
    chomp;
    @line = split(/\t/, $_);

    if ($line[1]==$unique_sites[$i]) {
	$i++;
	next;
    } else {
	print OUTPUT "$line[0]\t$line[1]\n";
    } #end else
    
} #end while

print STDERR "\tdone.\n";

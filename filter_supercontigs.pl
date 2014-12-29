#! /usr/bin/perl

$map = $ARGV[0];
$super = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide map file and super contig file\n\n";
    die;
} #end unless

open(MAP, $map);
open(SC, $super);
open(OUTPUT, ">$map.filtered");

@sc = ();
@sc_start = ();
@sc_end = ();
$block_header = (<SC>);

print STDERR "Reading in synteny blocks...";

while (<SC>) {
    chomp;
    @sc_array = split(/\s+/, $_);
    push(@sc, $sc_array[0]);
    push(@sc_start, $sc_array[3]/1000000);
    push(@sc_end, $sc_array[4]/1000000);
} #end while

$i=0;
$first_snp_between_contigs = 0;
$first_snp_in_new_contig = 0;

$header = (<MAP>);
chomp $header;
@map = ();
print OUTPUT "Synteny_block\t$header\tfilter\n";
print STDERR "done.\nReading in map file...";

while (<MAP>) {
    chomp;
    @line_item = split(/\t/,$_);
    push(@map, $_);
    push(@map_coord, $line_item[1]);
} #end while

print STDERR "done.\nProcessing map file for synteny blocks...";
$snp_counter = 0;
$num_snps_between = 0;

for ($k=0; $k<=$#map; $k++) {

    @input_array = split("\t", $map[$k]);

    if ($k==0) {
	$current_start = $map_coord[$k+50];
	for ($j=$k; $j<=$#map; $j++) {
	    if ($map_coord[$j]>=$sc_end[$i]) {
		$current_end = $map_coord[$j-50];
		last;
	    } #end if
	} #end for
#	print STDERR "Block: $sc[$i], start: $sc_start[$i]; 50+start: $current_start; end: $sc_end[$i]; 50-end: $current_end\n";
    } #end if

    if ($input_array[1]>=$sc_start[$i] && $input_array[1]<$current_start) { #within 50 SNPs of beginning of synteny block
	print OUTPUT "$sc[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t4\n";
	$snp_counter++;
    } elsif  ($input_array[1]>=$current_start && $input_array[1]<$current_end && $input_array[4]==100 && $input_array[5]==100 && $input_array[6]==100) { #before contig boundary end, filter1
	print OUTPUT "$sc[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t0\t0\t0\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t1\n";
	#print OUTPUT "$sc[$i]\t$map[$k]\t1\n";
	$snp_counter++;
    } elsif  ($input_array[1]>=$current_start && $input_array[1]<$current_end && $input_array[4]==150 && $input_array[5]==150 && $input_array[6]==150) { #before contig boundary end, filter2
	print OUTPUT "$sc[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t0\t0\t0\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t2\n";
	#print OUTPUT "$sc[$i]\t$map[$k]\t2\n";
	$snp_counter++;
    } elsif  ($input_array[1]>=$current_start && $input_array[1]<$current_end) { #before contig boundary end
	print OUTPUT "$sc[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t0\n";
	$snp_counter++;
    } elsif ($input_array[1]>=$current_end && $input_array[1]<=$sc_end[$i]) {#within 50 SNPs of end of synteny block
	print OUTPUT "$sc[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t4\n";
	$snp_counter++;
#    } elsif ($input_array[1]>=$sc_end[$i] && $input_array[1]<=$sc_start[$i+1] && $first_snp_between_contigs == 0) { #between contigs, first interval
#	$first_snp_between_contigs = 1;
#	$first_snp_in_new_contig = 1;
#	@first_interval_array = @input_array;	#need to gather coordinate info to collapse with last between contig
#	$num_snps_between++;
#    } elsif ($input_array[1]>=$sc_end[$i] && $input_array[1]<$sc_start[$i+1] && $first_snp_between_contigs == 1) { #between contigs, subsequent intervals
#	$num_snps_between++;
#	next;
    } elsif ($input_array[1]>=$sc_start[$i+1]) { #coordinates now within next super contig, first interval
	$i++;
	$current_start=$map_coord[$k+50];
	for ($j=$k; $j<=$#map; $j++) {
	    if ($map_coord[$j]>=$sc_end[$i]) {
		$current_end = $map_coord[$j-50];
		last;
	    } elsif ($sc[$i]==$sc[$#sc] && $map_coord[$j]>=$map_coord[$#map]) {
		$current_end = $map_coord[$j-50];
		last;
	    } #end if
	} #end for
#	print STDERR "num_snps: $snp_counter; num_between: $num_snps_between\nBlock: $sc[$i], start: $sc_start[$i]; 50+start: $current_start; end: $sc_end[$i]; 50-end: $current_end\n";
	$snp_counter=0;
	$num_snps_between=0;
	$first_snp_between_contigs = 0;
	$first_snp_in_new_contig = 0;
#	@last_interval_array = @input_array;
	#need to print collapsed interval
#	$new_midpoint = ($first_interval_array[1] + $last_interval_array[2])/2;
#	$hg18_new_midpoint = ($first_interval_array[8] + $last_interval_array[9])/2;
#	print OUTPUT "-\t$first_interval_array[0]\t$first_interval_array[1]\t$last_interval_array[2]\t$new_midpoint\t0\t0\t0\t$first_interval_array[7]\t$first_interval_array[8]\t$last_interval_array[9]\t$hg18_new_midpoint\t3\n";
#	@last_interval_array = ();
#	@first_interval_array = ();
    } #end elsif
} #end for

#print STDERR "num_snps: $snp_counter; num_between: $num_snps_between\n";
print STDERR "done.\n";

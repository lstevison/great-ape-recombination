#! /usr/bin/perl

#Calculates mean rate at predicted hotspots using LD map
#May 2, 2013
#Modified: July 30, 2014

$map_file = $ARGV[0];
$bed_file = $ARGV[1];
$output = $ARGV[2];
$chr_column_inmapfile = $ARGV[3];
$rate_column = $ARGV[4];
$map_format= $ARGV[5];

unless ($#ARGV==5 && $map_format=~/(bp|Mb)/) {
    print STDERR "Program reduces a map file based on a set of genomic coordinates in BED format.\nPlease provide on command line input map file, bed file, output filename, start column in the map file for the chromosome, column for the recombination rate, and the format of the map file (bp or Mb)\n\n";
    die;
} #end unless

open(MAP, $map_file);
open(BED, $bed_file);
open(OUTPUT, ">$output");

%map = (
"chr1" => undef,
"chr2" => undef,
"chr2a" => undef,
"chr2b" => undef,
"chr3" => undef,
"chr4" => undef,
"chr5" => undef,
"chr6" => undef,
"chr7" => undef,
"chr8" => undef,
"chr9" => undef,
"chr10" => undef,
"chr11" => undef,
"chr12" => undef,
"chr13" => undef,
"chr14" => undef,
"chr15" => undef,
"chr16" => undef,
"chr17" => undef,
"chr18" => undef,
"chr19" => undef,
"chr20" => undef,
"chr21" => undef,
"chr22" => undef);

@chr = ();
@start_pos = ();
@end_pos = ();

print STDERR "Now reading in bed input file...";

while (<BED>) {
    chomp;
    @input_array = split(/\s/, $_);
#    push(@chr, "chr" . $input_array[0]);
    push(@chr,$input_array[0]);
    push(@start_pos, $input_array[1]);
    push(@end_pos, $input_array[2]);
} #end while

print STDERR "First BED entry: $chr[0]: $start_pos[0]-$end_pos[0]\n";

print STDERR "done.\nNow reading in map file...";

while (<MAP>) {
    chomp;
    @map_array = split(/\s/, $_);

    if (defined($map_array[10]) && $map_array[12]!=0) {
	next;
    } #end if

    $map_chr = $map_array[$chr_column_inmapfile];
    if ($map_format=~/Mb/) {
#	if ($map_chr=~/chr2b/) {
#	    $map_start = ($map_array[$chr_column_inmapfile+1]-111.609048)*1000000;
#	    $map_end = ($map_array[$chr_column_inmapfile+2]-111.609048)*1000000;
#	} else {
	    $map_start = $map_array[$chr_column_inmapfile+1]*1000000;
	    $map_end = $map_array[$chr_column_inmapfile+2]*1000000;
#	} #end else
    } elsif ($map_format=~/bp/) {
	$map_start = $map_array[$chr_column_inmapfile+1];
	$map_end = $map_array[$chr_column_inmapfile+2];
    } #end else

    $rate = $map_array[$rate_column];
    $value = "$map_start:$map_end:$rate";
    push @{$map{$map_chr}}, $value;

} #end while

print STDERR "done.\n Now reducing map based on $#chr coordinates in bed file...\n";
print OUTPUT "chr\tstart\tend\tpeak_rho/kb\tavg_rho/kb\n";

for ($a=0; $a<=$#chr; $a++) {
#    print STDERR "Now on bed interval $a...\n";
    $total_rate = 0;
    $total_distance = 0;
    @peak_rate=();

    if ($chr[$a] ne $chr[$a-1] || $a==0) {
	@chr_matched_map = @{$map{$chr[$a]}};
#	@sorted_matched_map = sort @chr_matched_map;
	print STDERR "$chr[$a]\t$chr_matched_map[0], $chr_matched_map[1], $chr_matched_map[2], $chr_matched_map[3], $chr_matched_map[4], $chr_matched_map[5]\n";
	@map_start = ();
	@map_end = ();
	@rec_rate = ();

	for ($z=0; $z<=$#chr_matched_map; $z++) {
	    @data = split(":",$chr_matched_map[$z]);
	    push(@map_start, $data[0]);
	    push(@map_end, $data[1]);
	    push(@rec_rate, $data[2]);
	} #end for
    } #end if

    for ($b=0; $b<=$#map_start; $b++) {

#	$half_length=(($map_end[$b]-$map_start[$b])/2);
#	$int_min=$start_pos[$a]-$half_length;
#	$int_max=$end_pos[$a]-$half_length;

	if ($map_start[$b] >= $end_pos[$a]) {  #map coordinate after interval
	    last;
	} elsif ($map_end[$b] <= $start_pos[$a]) { #map coordinate before interval, keep going
	    next;
	} else {
#if (($map_start[$b] >= $start_pos[$a] && $map_end[$b] <= $end_pos[$a]) || ($map_start[$b] >= $start_pos[$a] && $int_max <= $end_pos[$a]) || ($int_min <= $start_pos[$a] && $map_end[$b] <= $end_pos[$a])) {  #map coordinate within interval
#	    print OUTPUT "$chr[$a]\t$map_start[$b]\t$map_end[$b]\t$rec_rate[$b]\t$block_number\n";
	    $distance = $map_end[$b]-$map_start[$b]+1;
	    $kb_distance = $distance/1000;
	    $total_rate += $rec_rate[$b]*$kb_distance;
	    $total_distance += $kb_distance;
	    push(@peak_rate, $rec_rate[$b]);
	    $max_b = $b;
	    $count+=1;
	} #end elsif
    } #end b-for
    if ($total_rate==0 || $total_distance==0) {
	$avg_rate = 0;
	print OUTPUT "$chr[$a]\t$start_pos[$a]\t$end_pos[$a]\tn/a\tn/a\tn/a\tn/a\n";
#	print STDERR "No matches at $chr[$a]:$start_pos[$a]-$end_pos[$a]\n";
    } else {
	@sorted_array = sort {$a <=> $b} @peak_rate;
	$avg_rate = $total_rate / $total_distance;
	print OUTPUT "$chr[$a]\t$start_pos[$a]\t$end_pos[$a]\t$sorted_array[$#sorted_array]\t$avg_rate\t$count\t$sorted_array[0]\n";
#	print STDERR "Intervals: $count; Max: $sorted_array[$#sorted_array]; Min: $sorted_array[0]; Average: $avg_rate\n";
    } #end

#    splice(@map_start,0,$max_b);
#    splice(@map_end,0,$max_b);
#    splice(@rec_rate,0,$max_b);
    $max_b=0;
    $count=0;
} #end a-for

print STDERR "done.\n";

#! /usr/bin/perl

#February 27, 2012
#Modified November 17, 2014

$bed = $ARGV[0];
$map = $ARGV[1];
$output = $ARGV[2];
$Ne = $ARGV[3];

unless ($#ARGV==3) {
    print STDERR "Please specify on command line bed input file, input map file, output filename, and corresponding Ne value.\n";
    die;
} #end unless

unless (-e $bed) {
    print STDERR "Error: $bed file does not exist!\n";
    die;
} #end unless

unless (-e $map) {
    print STDERR "Error: $map file does not exist!\n";
    die;
} #end if

open(BED, $bed);
open(INPUT, $map);
open(OUTPUT, ">$output");

print STDERR "Running hotspot calculations on map input $map and hotspots $bed\n";

%map = (
"chr1" => undef,
"chr2" => undef,
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
@hotspot_chromosome = ();
@positions = ();
@distance = ();
@rho = ();

while (<BED>) {
    chomp;
    @bed_array = split(/\t/, $_);
    push(@hotspot_chromosome, $bed_array[0]);
    $midpoint_position =  (($bed_array[2]+$bed_array[1]) / 2)/1000;
    push(@positions, $midpoint_position);
} #end while

print STDERR "Finished reading in $bed file\n";

$header=(<INPUT>);

while (<INPUT>) {
    chomp;
    @input_array = split(/\s+/, $_);

    $start_coordinate_kb = $input_array[1]/1000; #coordinates in bp
    $end_coordinate_kb = $input_array[2]/1000;
    $rec_rate = ($input_array[3] * 100000)/(4*$Ne); #rate in rho/kb, convert to cM/Mb
    $value = "$start_coordinate_kb:$end_coordinate_kb:$rec_rate";
    push @{$map{$input_array[0]}}, $value;

} #end while

print STDERR "Finished reading in $map file\n";

$number_hotspots = $#positions + 1;
print STDERR "Now performing calculations of rates near $number_hotspots known hotspots...\n";

for ($a=0; $a<=$#positions; $a++) {

    $start = $positions[$a] - 20;
    $end = $positions[$a] + 20;
#    if ($hotspot_chromosome[$a] ne $hotspot_chromosome[$a+1] || $a==$#positions) {
#	$last_hotspot_on_chr = 1;	  
#    } else {
#	$last_hotspot_on_chr = 0;
#    } #end else
#    print STDERR "A: $a; Calculating rates near $hotspot_chromosome[$a]:$positions[$a]\n";

    if ($hotspot_chromosome[$a] ne $hotspot_chromosome[$a-1] || $a==0) {
	@chr_matched_map = @{$map{$hotspot_chromosome[$a]}};
	print STDERR "$hotspot_chromosome[$a]\t$chr_matched_map[0], $chr_matched_map[1], $chr_matched_map[2], $chr_matched_map[3], $chr_matched_map[4], $chr_matched_map[5]\n";
	@start_pos = ();
	@end_pos = ();
	@rec_rate = ();

	for ($z=0; $z<=$#chr_matched_map; $z++) {
	    @data = split(":",$chr_matched_map[$z]);
	    push(@start_pos, $data[0]);
	    push(@end_pos, $data[1]);
	    push(@rec_rate, $data[2]);
	} #end for

    } #end if

    for ($i=0; $i<=$#start_pos; $i++) {
	
	if ($end_pos[$i]<$start) { #map coordinate before whole region, keep going through the map
	    next;
	} elsif ($start_pos[$i]>$end) { #map coordinate after whole region and no more hotspots on this chromosome - end it now!
	    last;
	} else {  #otherwise the coordinate should be in the hotspot!
	    $distance_to_motif = $positions[$a] - (($start_pos[$i]+$end_pos[$i])/2);
#	    print STDERR "I: $i; Length of chr: $#start_pos; Found a match at $hotspot_chromosome[$a]:$positions[$a]. Start: $start_pos[$i]; End: $end_pos[$i]; Rec rate: $rec_rate[$i]\n";
	    push(@distance, $distance_to_motif);				
	    push(@rho, $rec_rate[$i]);
	    next;
	} #end else

    } #end b-for

} #end a-for

print STDERR "Finished calculations of rates near all known hotspots\n";

%hash = ();
@sorted_rho = ();
@sorted_distance = ();

for ($x=0; $x<=$#rho; $x++) {
    unless ( abs($distance[$x])>20 ) {
	$hash{$distance[$x]} = $rho[$x];
    } #end unless
} #end for

print STDERR "Now printing output file\n";

foreach $value (sort { $a <=> $b } keys %hash) {
    push(@sorted_rho, $hash{$value});
    push(@sorted_distance, $value);
    print OUTPUT "$value\t$hash{$value}\n";
} #end foreach


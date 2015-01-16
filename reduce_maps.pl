#! /usr/bin/perl

#reduces map files based on an input bed file
#May 2, 2013
#Modified: May 6, 2013

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
#@bonobo_ori = ();
#@chimp_ori = ();
#@gorilla_ori = ();
print STDERR "Now reading in bed input file...";

while (<BED>) {
    chomp;
    @input_array = split(/\t/, $_);
    push(@chr, $input_array[0]);
    push(@start_pos, $input_array[1]);
    push(@end_pos, $input_array[2]);
#    push(@bonobo_ori, $input_array[3]);
#    push(@chimp_ori, $input_array[4]);
#    push(@gorilla_ori, $input_array[5]);
} #end while

print STDERR "done.\nNow reading in map file...";

while (<MAP>) {
    chomp;
    @map_array = split(/\t/, $_);

    if (defined($map_array[10]) && $map_array[12]!=0) {
	next;
    } #end if

    $map_chr = $map_array[$chr_column_inmapfile];
    if ($map_format=~/Mb/) {
	$map_start = $map_array[$chr_column_inmapfile+1]*1000000;
	$map_end = $map_array[$chr_column_inmapfile+2]*1000000;
    } elsif ($map_format=~/bp/) {
	$map_start = $map_array[$chr_column_inmapfile+1];
	$map_end = $map_array[$chr_column_inmapfile+2];
    } #end else

    $rate = $map_array[$rate_column];
    $value = "$map_start:$map_end:$rate";
    push @{$map{$map_chr}}, $value;

} #end while

print STDERR "done.\n Now reducing map based on $#chr coordinates in bed file...\n";
print OUTPUT "Hg18_chr\thg18_start\thg18_end\trho/kb\tsynteny_block\n";

for ($a=0; $a<=$#chr; $a++) {
    print STDERR "Now on interval $a...\n";
    $block_number = $a + 1;

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
	if ($map_start[$b] >= $start_pos[$a] && $map_end[$b] <= $end_pos[$a]) {  #map coordinate within interval
	    print OUTPUT "$chr[$a]\t$map_start[$b]\t$map_end[$b]\t$rec_rate[$b]\t$block_number\n";
	} #end else
    } #end b-for

} #end a-for

print STDERR "done.\n";

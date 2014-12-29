#! /usr/bin/perl

#Program reads in synteny blocks final output and map
#and re-prints map with hg18 coordinates at end
#Modified February 19, 2013

$synteny = $ARGV[0];
$input_map = $ARGV[1];
$output = $ARGV[2];

unless ($#ARGV==2) {
    print STDERR "Please provide on command line synteny block input, genetic map, and output filename\n\n";
    die;
} #end unless

open(SYNTENY, $synteny);
open(MAP, $input_map);
open(OUTPUT, ">$output");

@hg18_start = ();
@hg18_chr = ();
@gorgor_start = ();
@block = ();

#print header of output
print OUTPUT "Chr\tGorgor3-start\tGorgor3-end\tGorgor3-midpoint\trho/kb\tL95\tU95\thg18-chr\thg18-start\thg18-end\thg18-midpoint\n";

print STDERR "Now reading in synteny file...";

while (<SYNTENY>) {
    chomp;
    @line_item = split(/\s+/, $_);
    chomp($line_item[2]);
    push(@block, $line_item[0]);
    $gorilla = $line_item[2]/1;
    push(@gorgor_start, $gorilla);
    push(@hg18_chr, $line_item[3]);
    $human = $line_item[4]/1;
    push(@hg18_start, $human);
} #end while

#get rid of first 100 lines
for ($z=0; $z<98; $z++) {
    shift @hg18_start;
    shift @gorgor_start;
    shift @hg18_chr;
    shift @block;
} #end for

$first_map_line = 1;
$last_block = $block[0];

print STDERR "done.\nNow reading in map file and matching $#gorgor_start gorilla coordinates to $#hg18_start human coordinates...\nNow on block ";
#print STDERR "0:$gorgor_start[0]; 1:$gorgor_start[1]; 2:$gorgor_start[2]; 3:$gorgor_start[3]; 4:$gorgor_start[4]; 5:$gorgor_start[5]; Length: $#gorgor_start.";

while(<MAP>) {
    chomp;
    @map_line = split(/\t/, $_);
    chomp($map_line[1]);
    $gorgor_start_current = $map_line[1]*1000000;
    $gorgor_start_new = sprintf "%.0f", $gorgor_start_current;

    if ($gorgor_start_new==$gorgor_start_last) {
	print STDERR "Multiple map lines with same coordinates at $gorgor_start_current\n";
	next;
    } #end

    $gorgor_start_last = $gorgor_start_new;
    $chr = $map_line[0];
    $j=0;

    while ($j<=$#hg18_start) {
	if ($gorgor_start_new==$gorgor_start[$j] && $first_map_line==0) {
	    $hg18_start_new = $hg18_start[$j]/1000000;
	    $midpoint = ($last_hg18_start+ $hg18_start_new)/2;
	    print OUTPUT "$hg18_start_new\t$midpoint\n$chr\t$map_line[1]\t$map_line[2]\t$map_line[3]\t$map_line[4]\t$map_line[5]\t$map_line[6]\t$hg18_chr[$j]\t$hg18_start_new\t";
	    $last_hg18_start = $hg18_start_new;
	    $last_block = $block[$j];
	    splice (@hg18_start, $j, 1);
	    splice (@gorgor_start, $j, 1);
	    splice (@hg18_chr, $j, 1);
	    splice (@block, $j, 1);
#	    print STDERR "J: $j; size: $#hg18_chr ";
	    last;
	} elsif ($gorgor_start_new==$gorgor_start[$j] && $first_map_line==1) {
	    $hg18_start_new = $hg18_start[$j]/1000000;
	    print OUTPUT "$chr\t$map_line[1]\t$map_line[2]\t$map_line[3]\t$map_line[4]\t$map_line[5]\t$map_line[6]\t$hg18_chr[$j]\t$hg18_start_new\t";
	    print STDERR "$block[$j] ";
#	    print STDERR "Gorilla start: $gorgor_start[$j], Human start: $hg18_start[$j]...\n";
#	    print STDERR "J: $j; size: $#hg18_chr ";
	    $last_hg18_start = $hg18_start_new;
	    $last_block = $block[$j];
	    $first_map_line = 0;
	    splice (@hg18_start, $j, 1);
	    splice (@gorgor_start, $j, 1);
	    splice (@hg18_chr, $j, 1);
	    splice (@block, $j, 1);
	    last;
	} elsif ($last_block==$block[$j] && $first_map_line==0) { #non-match, but still in same block, could be MAF site
	    $maf = 0;
	    $edge = 0;
	    for ($i=$j; $i<=$#hg18_start; $i++) {
		if ($gorgor_start_new==$gorgor_start[$i] && $last_block==$block[$i]) { #test if matching site in current block, meaning it's a MAF site, otherwise edge sites
		    $maf = 1;
		    $num_maf = $i;
		    last;
		} elsif ($gorgor_start_new==$gorgor_start[$i] && $last_block!=$block[$i]) { #test if matching site not in current block, meaning it's an edge site
		    $edge = 1;
		    $num_sites_to_next = $i;
#		    print STDERR "At map site $gorgor_start_new found potential end of block: $gorgor_start[$j]; J: $j; Length array before: $#hg18_start (first element $gorgor_start[0]) ";
		    $hg18_start_new = $hg18_start[$j]/1000000;
		    $midpoint = ($last_hg18_start+ $hg18_start_new)/2;
		    print OUTPUT "$hg18_start_new\t$midpoint\n";
		    $first_map_line = 1;
		    last;
		} #end elsif
	    } #end for

	    $splice_counter = 0;

	    if ($maf==1) { #site found meaning it was removed for MAF cutoff, need to remove from list too
#		print STDERR "At map site $gorgor_start_new found maf site in block file: $gorgor_start[$j], J: $j; Num sites: $num_maf; Length array before: $#hg18_start (first element $gorgor_start[0]) ";
		while ($splice_counter<$num_maf) {
		    splice (@hg18_start, $j, 1);
		    splice (@gorgor_start, $j, 1);
		    splice (@hg18_chr, $j, 1);
		    splice (@block, $j, 1);	    
		    $splice_counter++;
		} #end while
#		print STDERR "After: $#hg18_start (first element $gorgor_start[0])\n";
	    } elsif ($edge==1) { #not a MAF site, so end of block, capture end site
		while ($splice_counter<$num_sites_to_next) {
		    splice (@hg18_start, $j, 1);
		    splice (@gorgor_start, $j, 1);
		    splice (@hg18_chr, $j, 1);
		    splice (@block, $j, 1);	    
		    $splice_counter++;
		} #end while
#		print STDERR "After: $#hg18_start (first element $gorgor_start[0])\n";
		print STDERR "done\nNow on block ";
#	    } else {
#		print STDERR "Map file site: $gorgor_start_new / Synteny file site: $gorgor_start[$j] is not an maf OR edge site, but does not appear to have a match in synteny file\n\n";
#		die;
	    } #end else
	} else {
	    $j++;
#	    print STDERR "Match not found for map site $gorgor_start_new; J: $j; Gor: $gorgor_start[$j]; Human: $hg18_start[$j]; Last hg18: $last_hg18_start\n";
	} #end elsif
    } #end while
} #end while

#print STDERR "Size of array at end: $#hg18_start";
$hg18_start_new = $hg18_start[$j]/1000000;
$midpoint = ($last_hg18_start + $hg18_start_new)/2;
print OUTPUT "$hg18_start_new\t$midpoint\n";

print STDERR "done.\n";

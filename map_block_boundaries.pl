#! /usr/bin/perl

#program reads through each chromosome, and records the start, stop and orientation for each syntenic block
#Created: July 16, 2013
#Modified: July 24, 2014

$chr_list = 'chr_list.txt';
$dir=$ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide input directory and output BED filename on command line\n\n";
    die;
} #end unless

unless (-e $chr_list) {
    print STDERR "Error: File $chr_list cannot be found.\n\n";
    die;
} #end 

open(LIST, $chr_list);
@chr_list = ();

while(<LIST>) {
    chomp;
    push(@chr_list, $_);
} #end while

open(OUTPUT, ">$output");

print OUTPUT "Non-human chr\tstart\tstop\tHuman chr\tstart\tstop\tOrientation\tBlock\tOld block\tNum sites\n";

$block=1;

for ($i=0; $i<=$#chr_list; $i++) {

    $map = $dir . $chr_list[$i] . '-map.hg18.txt.filtered';

    if (-e $map) {
	print STDERR "Reading through $chr_list[$i] map file and recording coordinates...\n";
    } else {
	print STDERR "Error: file $map does not exist!\n\n";
	die;
    } #end 


    open(MAP, $map);
    $header=(<MAP>);
    $new_block = 1;
    $block2=1;
    $count_sites = 0;
    
    while(<MAP>) {
	chomp;
	@input_line = split(/\t/, $_);
	
	if ($input_line[0]==$block2 && $new_block==1) {
	    $n_start = $input_line[2]*1000000;
	    $n_end = $input_line[3]*1000000;
	    $n_chr = $input_line[1];
	    $h_start = $input_line[9]*1000000;
	    $h_end = $input_line[10]*1000000;
	    $h_chr = $input_line[8];
	    $new_block=0;
	    $count_sites=1;
	    next;
	} elsif ($input_line[0]!=$block2 && $new_block==0 && $count_sites>1) {
	    if ($last_ori=~/plus/) {
		$h_min=$h_start;
		$h_max=$last_h_end;
		$n_min=$n_start;
		$n_max=$last_n_end;
	    } elsif ($last_ori=~/minus/) {
		$n_min=$n_start;
		$n_max=$last_n_end;
		$h_min=$last_h_end;
		$h_max=$h_start;
	    } #end if

	    print OUTPUT "$n_chr\t$n_min\t$n_max\t$h_chr\t$h_min\t$h_max\t$last_ori\t$block\t$block2\t$count_sites\n";
	    $block++;
	    $count_sites = 1;
	    $block2 = $input_line[0];
	    $n_start = $input_line[2]*1000000;
	    $n_end = $input_line[3]*1000000;
	    $n_chr = $input_line[1];
	    $h_start = $input_line[9]*1000000;
	    $h_end = $input_line[10]*1000000;
	    $h_chr = $input_line[8];
	    next;
	} elsif ($input_line[0]==($block2+1)) {
	    print STDERR "Discarding block $block2, which has only $count_sites site.\n";
	    $block2 = $input_line[0];
	    $n_start = $input_line[2]*1000000;
	    $n_end = $input_line[3]*1000000;
	    $n_chr = $input_line[1];
	    $h_start = $input_line[9]*1000000;
	    $h_end = $input_line[10]*1000000;
	    $h_chr = $input_line[8];
	    $new_block=0;
	    next;
	} elsif ($new_block==0) {
	    $count_sites++;
	    $last_h_start = $input_line[9]*1000000;
	    $last_h_end = $input_line[10]*1000000;
	    $last_n_start = $input_line[2]*1000000;
	    $last_n_end = $input_line[3]*1000000;
	    $last_ori = $input_line[13];
	    next;
	} #end if

	print STDERR "Error at line: $_\n\n";
	$block2 = $input_line[0];

    } #end while 

    if ($last_ori=~/plus/) {
	$h_min=$h_start;
	$h_max=$last_h_end;
	$n_min=$n_start;
	$n_max=$last_n_end;
    } elsif ($last_ori=~/minus/) {
	$n_min=$n_start;
	$n_max=$last_n_end;
	$h_min=$last_h_end;
	$h_max=$h_start;
    } #end if
    $count_sites++;
    print OUTPUT "$n_chr\t$n_min\t$n_max\t$h_chr\t$h_min\t$h_max\t$last_ori\t$block\t$block2\t$count_sites\n";
    $block++;

} #end for

print STDERR "done.\n";

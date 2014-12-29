#! /usr/bin/perl

#July 24, 2012
#Modified October 25, 2012
#Program reads in output of 'convert-cood-gorgor2hg18.pl'
#defines synteny blocks

print STDERR "Running program to define syteny blocks\n";

$input = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide on command line input and output filenames\n\n";
    die;
} #end unless

open(INPUT, $input);

@gor_chr = ();
@gor_coord = ();
@hg18_chr = ();
@hg18_coord = ();

print STDERR "Reading input...";

while(<INPUT>) {
    chomp;
    @input_line = split(/\t/, $_);
    push(@gor_chr, $input_line[0]);
    push(@gor_coord, $input_line[1]);
    push(@hg18_chr, $input_line[2]);
    push(@hg18_coord, $input_line[3]);
} #end while

print STDERR "done.\n";

@synteny_block = ();
$region = 1;
$first_line_block = 1;
$output2 = $output . "2";
$distance_cutoff=0;

open(OUTPUT, ">$output");
open(OUTPUT2, ">$output2");
print OUTPUT "Region\tGorilla_chr\tGorilla_coord\tHuman_chr\tHuman_coord\tOrientation\n"; #print header
print OUTPUT2 "Region\tGorilla_chr\tGorilla_start\tGorilla_end\thg18_chr\thg18_start\thg18_end\tOrientation\t\#SNPs\t\#Bad_snps\tDistance break?\n";
print STDERR "Now processing synteny blocks...";

LOOP:

$bad_snp_counter = 0;
$snp_counter = 0;
$big_jump=0;

for ($i=0; $i<=$#gor_coord; $i++) {

    $hg18_distance = abs($hg18_coord[$i+1] - $hg18_coord[$i]);
    $gor_distance = abs($gor_coord[$i+1] - $gor_coord[$i]);

    if ($gor_distance>=1000000 || $hg18_distance>=1000000) {
	$big_jump=1;
    } #end if

    if ($gor_coord[$i]< $gor_coord[$i+1] && $hg18_coord[$i]< $hg18_coord[$i+1] && $hg18_chr[$i] eq $hg18_chr[$i+1] && $gor_chr[$i] eq $gor_chr[$i+1] && $big_jump==0) {
	$orientation = "plus";
    } elsif ($gor_coord[$i]< $gor_coord[$i+1] && $hg18_coord[$i]> $hg18_coord[$i+1] && $hg18_chr[$i] eq $hg18_chr[$i+1] && $gor_chr[$i] eq $gor_chr[$i+1] && $big_jump==0) {
	$orientation = "minus1";
    } elsif ($gor_coord[$i]> $gor_coord[$i+1] && $hg18_coord[$i]< $hg18_coord[$i+1] && $hg18_chr[$i] eq $hg18_chr[$i+1] && $gor_chr[$i] eq $gor_chr[$i+1] && $big_jump==0) {
	$orientation = "minus2";
    } else {
	$orientation = "bad";
    } #end else

    if ($i==0 && $big_jump==0) {
	$last_orientation = $orientation;
	$first_block=1;
	next;
    } elsif ($i==0 && $big_jump==1) {
	shift(@gor_coord);
	shift(@gor_chr);
	shift(@hg18_coord);
	shift(@hg18_chr);
	goto LOOP;
    } #end if

    if ($last_orientation eq $orientation && $first_line_block==1 && $gor_distance<=50000 && $hg18_distance<=50000) {
	$first_line_block = 0;
	$snp_counter=0;

	if ($bad_snp_counter >=1 && $last_orientation=~/plus/ && $gor_coord[$i]> $gor_coord[$i-1] && $hg18_coord[$i]> $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1] && $distance_cutoff==0 && $old_hg18_distance<=50000 && $old_gor_distance<=50000) {
	    $first_synteny_line = "$region\t$gor_chr[$i-1]\t$gor_coord[$i-1]\t$hg18_chr[$i-1]\t$hg18_coord[$i-1]\t$orientation";
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tplus";
	    push(@synteny_block, $first_synteny_line, $synteny_line);
	    $snp_counter+=2;
	    if ($region>=2 && $bad_snp_counter>1) {
		for ($h=$bad_snp_counter;$h>1; $h--) {
		    print OUTPUT "-\t$gor_chr[$i-$h]\t$gor_coord[$i-$h]\t$hg18_chr[$i-$h]\t$hg18_coord[$i-$h]\t$orientation\n";
#		    print STDERR "Made it into loop; Region: $region; H: $h; Bad snps: $bad_snp_counter; $gor_coord[$i-$h]\n";
		} #endfor
	    } #end if
	} elsif ($first_block==1) {
	    $first_synteny_line = "$region\t$gor_chr[$i-1]\t$gor_coord[$i-1]\t$hg18_chr[$i-1]\t$hg18_coord[$i-1]\t$orientation";
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\t$orientation";
	    push(@synteny_block, $first_synteny_line, $synteny_line);
	    $snp_counter+=2;
	    $first_block=0;
	} elsif ($bad_snp_counter >=1 && $last_orientation=~/minus1/ && $gor_coord[$i]> $gor_coord[$i-1] && $hg18_coord[$i]< $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1] && $distance_cutoff==0 && $old_hg18_distance<=50000 && $old_gor_distance<=50000) {
	    $first_synteny_line = "$region\t$gor_chr[$i-1]\t$gor_coord[$i-1]\t$hg18_chr[$i-1]\t$hg18_coord[$i-1]\t$orientation";
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\t$orientation";
	    push(@synteny_block, $first_synteny_line, $synteny_line);
	    $snp_counter+=2;
	    if ($region>=2 && $bad_snp_counter>1) {
		for ($h=$bad_snp_counter;$h>1; $h--) {
		    print OUTPUT "-\t$gor_chr[$i-$h]\t$gor_coord[$i-$h]\t$hg18_chr[$i-$h]\t$hg18_coord[$i-$h]\t$orientation\n";
#		    print STDERR "Made it into loop; Region: $region; H: $h; Bad snps: $bad_snp_counter; $gor_coord[$i-$h]\n";
		} #endfor
	    } #end if
	} elsif ($bad_snp_counter >=1 && $last_orientation=~/minus2/ && $gor_coord[$i]< $gor_coord[$i-1] && $hg18_coord[$i]> $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1] && $distance_cutoff==0 && $old_hg18_distance<=50000 && $old_gor_distance<=50000) {
	    $first_synteny_line = "$region\t$gor_chr[$i-1]\t$gor_coord[$i-1]\t$hg18_chr[$i-1]\t$hg18_coord[$i-1]\t$orientation";
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\t$orientation";
	    push(@synteny_block, $first_synteny_line, $synteny_line);
	    $snp_counter+=2;
	    if ($region>=2 && $bad_snp_counter>1) {
		for ($h=$bad_snp_counter;$h>1; $h--) {
		    print OUTPUT "-\t$gor_chr[$i-$h]\t$gor_coord[$i-$h]\t$hg18_chr[$i-$h]\t$hg18_coord[$i-$h]\t$orientation\n";
#		    print STDERR "Made it into loop; Region: $region; H: $h; Bad snps: $bad_snp_counter; $gor_coord[$i-$h]\n";
		} #endfor
	    } #end if
	} else {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\t$orientation";
	    push(@synteny_block, $synteny_line);
	    $snp_counter++;
	    if ($region>=2 && $bad_snp_counter>=1) {
		for ($h=$bad_snp_counter;$h>=1; $h--) {
		    print OUTPUT "-\t$gor_chr[$i-$h]\t$gor_coord[$i-$h]\t$hg18_chr[$i-$h]\t$hg18_coord[$i-$h]\t$orientation\n";
#		    print STDERR "Made it into loop; Region: $region; H: $h; Bad snps: $bad_snp_counter; $gor_coord[$i-$h]\n";
		} #endfor
	    } #end if
	} #end if
	$distance_cutoff=0;
    } elsif ($last_orientation eq $orientation && $first_line_block==0 && $gor_distance<=50000 && $hg18_distance<=50000) {
	$synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\t$orientation";
	push(@synteny_block, $synteny_line);
	$snp_counter++;
    } elsif ($i==$#gor_coord && $#synteny_block>=0) {
	if ($last_orientation=~/plus/ && $gor_coord[$i]> $gor_coord[$i-1] && $hg18_coord[$i]> $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1]) {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tplus";
	    push(@synteny_block, $synteny_line);
	} elsif ($last_orientation=~/minus1/ && $gor_coord[$i]> $gor_coord[$i-1] && $hg18_coord[$i]< $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1]) {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tminus1";
	    push(@synteny_block, $synteny_line);
	} elsif ($last_orientation=~/minus2/ && $gor_coord[$i]< $gor_coord[$i-1] && $hg18_coord[$i]> $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1]) {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tminus2";
	    push(@synteny_block, $synteny_line);
	} #end if

	for ($j=0;$j<=$#synteny_block; $j++) {
	    print OUTPUT "$synteny_block[$j]\n";
	} #endfor

	@first_line = split(/\t/, $synteny_block[0]);
	@last_line = split(/\t/, $synteny_block[$#synteny_block]);
	print OUTPUT2 "$first_line[0]\t$first_line[1]\t$first_line[2]\t$last_line[2]\t$first_line[3]\t$first_line[4]\t$last_line[4]\t$first_line[5]\t$snp_counter\t$bad_snp_counter\t$distance_cutoff\n";
    } elsif ($#synteny_block>=0) {

	$bad_snp_counter=0; 
	$old_hg18_distance = $hg18_distance;
	$old_gor_distance = $gor_distance;

	if ($big_jump==0 && $last_orientation eq $orientation) {
	    if ($gor_distance>=50000 || $hg18_distance>=50000) {
		$distance_cutoff=1;
	    } #end if
	} #end if

	if ($last_orientation=~/plus/ && $gor_coord[$i]> $gor_coord[$i-1] && $hg18_coord[$i]> $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1] && $gor_coord[$i]<$gor_coord[$i+1] && $hg18_coord[$i]<$hg18_coord[$i+1] && $big_jump==0) {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tplus";
	    push(@synteny_block, $synteny_line);
	} elsif ($last_orientation=~/minus1/ && $gor_coord[$i]> $gor_coord[$i-1] && $hg18_coord[$i]< $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1] && $gor_coord[$i]<$gor_coord[$i+1] && $hg18_coord[$i]>$hg18_coord[$i+1] && $big_jump==0) {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tminus1";
	    push(@synteny_block, $synteny_line);
	} elsif ($last_orientation=~/minus2/ && $gor_coord[$i]< $gor_coord[$i-1] && $hg18_coord[$i]> $hg18_coord[$i-1] && $hg18_chr[$i] eq $hg18_chr[$i-1] && $gor_chr[$i] eq $gor_chr[$i-1] && $gor_coord[$i]>$gor_coord[$i+1] && $hg18_coord[$i]<$hg18_coord[$i+1] && $big_jump==0) {
	    $synteny_line = "$region\t$gor_chr[$i]\t$gor_coord[$i]\t$hg18_chr[$i]\t$hg18_coord[$i]\tminus2";
	    push(@synteny_block, $synteny_line);
	} else {
	    $bad_snp_counter++;
	} #end if

	for ($j=0;$j<=$#synteny_block; $j++) {
	    print OUTPUT "$synteny_block[$j]\n";
	} #endfor
	@first_line = split(/\t/, $synteny_block[0]);
	@last_line = split(/\t/, $synteny_block[$#synteny_block]);
	print OUTPUT2 "$first_line[0]\t$first_line[1]\t$first_line[2]\t$last_line[2]\t$first_line[3]\t$first_line[4]\t$last_line[4]\t$first_line[5]\t$snp_counter\t$bad_snp_counter\t$distance_cutoff\n";
	$region++;
	$snp_counter=0;
	@synteny_block = ();
	$first_line_block = 1;
	$big_jump=0;
#    } elsif ($#synteny_block==0) {
#	@first_line = split(/\t/, $synteny_block[0]);
#	print OUTPUT "-\t$first_line[1]\t$first_line[2]\t$first_line[3]\t$first_line[4]\t$first_line[5]\n";
#	$snp_counter=0;
#	@synteny_block = ();
#	$first_line_block = 1;
#	$big_jump=0;
#	$bad_snp_counter++;
    } else {
	$bad_snp_counter++; 
	$big_jump=0;
    } #end else

    $last_orientation = $orientation;

} #end for

print STDERR "done.\n";

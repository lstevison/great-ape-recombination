#! /usr/bin/perl

#Program to sort through 1kb regions within 100kb of hotspots and identify best matched cold spots
#Created: December 19, 2013
#Last modified: October 9, 2014

$hotspots_fasta=$ARGV[0];
$hotspots_withrates=$ARGV[1];
$input=$ARGV[2];
$output=$ARGV[3];
$genome_avg=$ARGV[4];
$peak_cutoff=$ARGV[5];
$width=$ARGV[6];

unless ($#ARGV==6) {
    print STDERR "Please provide on command line filename of fasta formatted hotspots, filename for hotspots with rates, a list of 1kb regions for the whole genome, an output filename, the genome average rate, the peak cutoff used for hotspots and the width to search for matching coldspots.\n\n";
    die;
} #end unless

open(HOTSPOTS, $hotspots_fasta);
open(RATES, $hotspots_withrates);
open(INPUT, $input);
$header=(<INPUT>);
open(OUTPUT, ">$output");

@hotspot_names = ();
@hotspot_sequences = ();
@hotspot_peaks=();
@hotspot_avg=();

print STDERR "Reading in hotspots...";

while (<HOTSPOTS>) {
    chomp;

    if (/>/) { 
        $seq_name = $'; #'
            push(@hotspot_names, $seq_name);
        if ($newseq==1) {push(@hotspot_sequences, $sequence);} 
        $sequence = "";
        $newseq=0;
        next; 
    } else {
        $sequence .= $_;
        $newseq=1; 
    } #end else    
} #end while
push(@hotspot_sequences, $sequence);
$h=0;
$header=(<RATES>);
while (<RATES>) {
    chomp;

    @rates=split(/\s+/,$_);
    $new_hotspot="$rates[0].$rates[1]:$rates[2]";

    if ($hotspot_names[$h]=~$new_hotspot) {
	push(@hotspot_peaks,$rates[3]);
	push(@hotspot_avg,$rates[4]);
	$h++;
    } #end if

} #end while

print STDERR "done reading in $#hotspot_names hotspots and $#hotspot_peaks rates.\nNow reading in coldspots for matching...";

@coldspot_matches_array=();
$l=0;

while (<INPUT>) {
    chomp;
    push(@coldspot_matches_array,$_);
} #end while

print STDERR "done reading in $#coldspot_matches_array coldspots for matching.\nNow matching hotspots to coldspots...";
print OUTPUT "Hotspot_chr\tHotspot_start\tHotspot_end\tHotspot_GC\tHotspot_N\tHotspot_avgRate\tHotspot_PeakRate\tColdspot_chr\tColdspot_start\tColdspot_end\tColdspot_GC\tColdspot_N\tColdspot_avgRate\tColdspot_PeakRate\n";

for ($i=0; $i<=$#hotspot_names; $i++) {

    $hotspot_rate=0;
    @input_array = split(/\.|:/, $hotspot_names[$i]);
#    print STDERR "$input_array[0]\t$input_array[1]\t$input_array[2]\n";
    if ($i==0) {
	$last_chr=$input_array[0];
    } #end if
    $hot_start=$input_array[1]+1;
    $hotspot_center = $input_array[1] + (($input_array[2]-$input_array[1])/2);
    $hotspot_length=sprintf("%.0f",($input_array[2]-$input_array[1])/1000);
    $new_width=(int($width/$hotspot_length))*$hotspot_length;

#    if ($i<=10) {print STDERR "Hotspot length: $hotspot_length, Width: $width; New width: $new_width\n";}
    if ($input_array[1] < $new_width) {    
	$max = $hotspot_center + $new_width + 1;
	$min = 0;
    } else {
	$max = $hotspot_center + $new_width + 1;
	$min = $hotspot_center - $new_width - 1;
    } #end else

    $a=($hotspot_sequences[$i]=~tr/Aa//);
    $t=($hotspot_sequences[$i]=~tr/Tt//);
    $g=($hotspot_sequences[$i]=~tr/Gg//);
    $c=($hotspot_sequences[$i]=~tr/Cc//);
    $total=$a+$t+$g+$c;
    $interval_size = $input_array[2]-$input_array[1]+1;
    $percent_N=sprintf("%.2f",(($interval_size-$total)/$interval_size)*100);
    if ($total != 0) {
	$gc= sprintf("%.2f", (($g+$c)/$total)*100 );
    } else {
	$gc=0;
    } #end else

    print OUTPUT "$input_array[0]\t$input_array[1]\t$input_array[2]\t$gc\t$percent_N\t$hotspot_avg[$i]\t$hotspot_peaks[$i]\t";
    @matched_coldspots = ();
    $previous_difference = 2;
    $first_match=0;

    for ($k=$l; $k<=$#coldspot_matches_array; $k+=$hotspot_length) {
	$highest_k=sprintf("%.0f", $k+$hotspot_length-1);
	$rate_count=0;
	$count=0;
	$coldpeak=0;
	$coldavg=0;
	$coldgc=0;
	$coldN=0;
	$start="";
	$end="";
	$cold_chr="";

	for ($b=$k; $b<=$highest_k; $b++) {
	    @region_array = split(/\s+/, $coldspot_matches_array[$b]);

	    if ($region_array[0]=~/^$input_array[0]$/) { 
		$cold_chr=$region_array[0];
		if ($b==$k) { $start=$region_array[1]; } #end if 
		if ($b==$highest_k) { $end=$region_array[2]; } #end if 

		if ($region_array[3]!~/n\/a/) {
		    $coldpeak+=$region_array[3];
		    $coldavg+=$region_array[4];
		    $rate_count++;
		} #end if
		$coldgc+=$region_array[7];
		$coldN+=$region_array[8];
		$count++;
	    } #end if
	} #end for

	unless (defined($end)) {
	    next;
	} #end unless

	unless (defined($cold_chr)) {
	    next;
	} #end unless

	if ($rate_count>0) { 
	    $coldspot_peak=$coldpeak/$rate_count;
	    $coldspot_avg=$coldavg/$rate_count;
	} else {
	    $coldspot_peak="n/a";
	    $coldspot_avg="n/a";
	} #end else

	if ($count>0) {
	    $coldspot_gc=$coldgc/$count;
	    $coldspot_N=$coldN/$count;
	} #end if

	$coldspot_center = $start + (($end-$start)/2);
#	print STDERR "$cold_chr\t$start\t$end\t$coldspot_peak\t$coldspot_avg\t$coldspot_gc\t$coldspot_N\n";
	if ($first_match==1) {
	    $new_l=$k-1;
	    $first_match=2;
	} #end 

	if ($input_array[0]!~/^$last_chr$/ && $cold_chr=~/^$last_chr$/) { #coldspot on last chr, but hotspot on new chr
#	    print STDERR "Hotspot chr: $input_array[0]; Last chr: $last_chr; Coldspot chr: $cold_chr\n";
	    next;
	} elsif ($coldspot_center >=$min && $coldspot_center <= $max && $cold_chr=~/^$input_array[0]$/) { #cold spot center within 100kb of hotspot

	    $coldspot_length=sprintf("%.0f", ($end-$start)/1000);
	    
	    unless ($coldspot_length==$hotspot_length) {
#		print STDERR "distance disparity, hotspot: $hotspot_length; coldspot: $coldspot_length\t$cold_chr:$start-$end\t$input_array[0]:$input_array[1]-$input_array[2]\n";
		next;
	    } #end if

	    if ($first_match==0) {
		$first_match=1;
	    } #end if

	    if ($coldspot_peak!~/n\/a/ && $coldspot_avg <= $genome_avg && $coldspot_peak <= $peak_cutoff && $$start!=$input_array[1] && $coldspot_N<=25) { #cold spot avg rate <=0.827 rho/kb (genome average), peak rate <=4.113 rho/kb (hotspot cutoff), not at hotspot, and percent N <=XX%
		$difference = abs($gc - $coldspot_gc);
#		print STDERR "Difference: $difference\n";
		if ($previous_difference > $difference) {
		    @matched_coldspots = ();
		    push(@matched_coldspots, "$cold_chr\t$start\t$end\t$coldspot_peak\t$coldspot_avg\t$coldspot_gc\t$coldspot_N");
		    $previous_difference = $difference;
		    next;
		} elsif ($previous_difference==$difference) {
		    push(@matched_coldspots, "$cold_chr\t$start\t$end\t$coldspot_peak\t$coldspot_avg\t$coldspot_gc\t$coldspot_N");
		} else {
		    next;
		} #end else
		
#	    } elsif ($start==$hot_start) { #print percentN, avg and peak rate at hotspot
#		print OUTPUT "$coldspot_N\t$coldspot_avg\t$coldspot_peak\t";
#		$hotspot_rate = 1;
	    } else {
		next;
	    } #end else
	} elsif ($end > $max) {
#	    if ($i<100) {print STDERR "L: $l; K:$k; Hotspot center: $hotspot_center; Min: $min; Max: $max; Coldspot center: $coldspot_center; Coldspot end: $region_array[2]; Matches: $#matched_coldspots\n";}
	    $l=$new_l;
	    last;
	} #end else
    } #end k-for

#    if ($hotspot_rate==0) {
#	print OUTPUT "n/a\tn/a\tn/a\t"; #if no rate for hotspot available, print 3 extra tabs so columns line up
#   } #end if

    if ($#matched_coldspots >= 1) {
	$previous_distance = $new_width;
	for ($z=0; $z<=$#matched_coldspots; $z++) {
	    @coldspot_array = split(/\s+/, $matched_coldspots[$z]);
	    $coldspot_center = $coldspot_array[1] + (($coldspot_array[2]-$coldspot_array[1])/2);
	    $distance = abs($hotspot_center - $coldspot_center);

	    if ($previous_distance > $distance) {
		@new_matched_coldspot = ();
		push(@new_matched_coldspot, $matched_coldspots[$z]);
		$previous_distance = $distance;
#		print STDERR "Distance to hotspot: $distance\n";
	    } #end if

	} #end for
	@coldspot_array = split(/\s+/, $new_matched_coldspot[0]);
	print OUTPUT "$coldspot_array[0]\t$coldspot_array[1]\t$coldspot_array[2]\t$coldspot_array[5]\t$coldspot_array[6]\t$coldspot_array[4]\t$coldspot_array[3]\tMultipleMatchesFound:$#matched_coldspots\n";
#	print STDERR "Multiple Matches found!\n";
    } elsif (defined($matched_coldspots[0])) {
	@coldspot_array = split(/\s+/, $matched_coldspots[0]);
	print OUTPUT "$coldspot_array[0]\t$coldspot_array[1]\t$coldspot_array[2]\t$coldspot_array[5]\t$coldspot_array[6]\t$coldspot_array[4]\t$coldspot_array[3]\n";
#	print STDERR "Match found!\n";
    } else {
#	print STDERR "no matches for @input_array\n";
	print OUTPUT "\t\t\t\t\t\t\tNoMatchesFound, $l\n";
    } #end if
    $last_chr=$input_array[0];
} #end for 

print STDERR "done.\n";

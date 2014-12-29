#! /usr/bin/perl

#Modified February 8, 2013

$map = $ARGV[0];
$synteny_blocks = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide map file and synteny block final2 file\n\n";
    die;
} #end unless

open(MAP, $map);
open(SB, $synteny_blocks);
open(OUTPUT, ">$map.filtered");

@sb = ();
@sb_start = ();
@sb_end = ();
@sb_ori = ();
$block_header = (<SB>);

print STDERR "Reading in synteny blocks...";

while (<SB>) {
    chomp;
    @sb_array = split(/\s+/, $_);
    push(@sb, $sb_array[0]);
    push(@sb_start, $sb_array[3]/1000000);
    push(@sb_end, $sb_array[4]/1000000);
    push(@sb_ori, $sb_array[10]);
} #end while

$i=0;

$header = (<MAP>);
chomp $header;
print OUTPUT "Synteny_block\t$header\tfilter\torientation\n";
print STDERR "done.\nReading in map file and processing for synteny blocks...";

while (<MAP>) {
    chomp;
    @input_array = split(/\t/,$_);

    if  ($input_array[1]>=$sb_start[$i] && $input_array[1]<=$sb_end[$i] && $input_array[4]==100 && $input_array[5]==100 && $input_array[6]==100) { #filter1 = high rho
	print OUTPUT "$sb[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\tna\tna\tna\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t1\t$sb_ori[$i]\n";
    } elsif  ($input_array[1]>=$sb_start[$i] && $input_array[1]<=$sb_end[$i] && $input_array[4]==150 && $input_array[5]==150 && $input_array[6]==150) { #filter2 = 50+/- high rho site
	print OUTPUT "$sb[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\tna\tna\tna\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t2\t$sb_ori[$i]\n";
    } elsif  ($input_array[1]>=$sb_start[$i] && $input_array[1]<=$sb_end[$i]) { #within synteny block i
	print OUTPUT "$sb[$i]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t0\t$sb_ori[$i]\n";
    } elsif ($input_array[1]>=$sb_start[$i+1]) { #coordinates now within next synteny block i+1
	print OUTPUT "$sb[$i+1]\t$input_array[0]\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[6]\t$input_array[7]\t$input_array[8]\t$input_array[9]\t$input_array[10]\t0\t$sb_ori[$i]\n";
	$i++;
    } #end elsif

} #end while

print STDERR "done.\n";

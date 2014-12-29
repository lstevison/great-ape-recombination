#! /usr/bin/perl

#May 8, 2012
#Modified August 16, 2012

#program reads in a position file and extracts coordinates
#in bed format
#of a specific number of SNPs in a window

$positions = $ARGV[0];	
$block = $ARGV[1];
$chr = $ARGV[2];

unless ($#ARGV==2) {
    print "Please provide filename of positions file, synteny block number and chromosome.\n\n";
    die;
} #end unless

open(INPUT, $positions);
open(OUTPUT, ">$chr.synteny_block.$block-input.BED");

$SNP_counter = 0;
@coordinates = ();
$first_SNP = 1;
$second_SNP = 0;

$header=(<INPUT>);

print STDERR "Processing synteny block: $block...\n";

while (<INPUT>) {
    chomp;
    $input_line = $_;
    @variant = split(/\s+/, $input_line);
    $SNP_counter++;
	
    if ($SNP_counter == 1 && $first_SNP == 1) {	 #prints first line of BED file
	print OUTPUT "$variant[0]\t" ;
	$first_SNP = 0;
	$first_interval = 1;
	$after_first_interval = 0;
    } elsif ($SNP_counter == 199 && $made_past_200 ==1 && $after_first_interval==0) {  #prints end position and next start position
	print OUTPUT "$variant[0]\n$start_next\t" ;
	$made_past_200 = 0;
	$second_SNP = 1;
	$after_first_interval = 1;
    } elsif ($SNP_counter == 199 && $made_past_200 ==1 && $second_SNP==1) {   	     #prints end position and next start position
	print OUTPUT "$variant[0]\n$start_next\t" ;
	$made_past_200 = 0;
	$second_SNP = 0;
    } elsif ($SNP_counter == 199 && $made_past_200 ==1 && $second_SNP==0) {	       #prints end position and next start position
	print OUTPUT "$variant[0]\n$start_next\t" ;
	$made_past_200 = 0;
    } elsif ($SNP_counter == 3800 && $after_first_interval == 1 && $made_it_past_200==0) {	     #holds next beginning
	$start_next = $variant[0];
	$SNP_counter=0;
	$made_past_200 =1;
    } elsif ($SNP_counter == 3801 && $first_interval == 1) {	      #holds next beginning
	$start_next = $variant[0];
	$SNP_counter=0;
	$made_past_200 =1;
	$first_interval = 0;
    } #end if
	
    $last_position = $variant[0];
    
} #end while

print OUTPUT "$last_position\n";

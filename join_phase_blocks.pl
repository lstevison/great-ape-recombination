#! /usr/bin/perl

#Program reads in PHASE outputs and pastes larger synteny blocks back together
#Creator: Laurie Stevison
#Last Modified: June 10, 2014

$chr=$ARGV[0];

unless ($#ARGV==0) {
    print STDERR "Please provide chromsome on command line\nEx. chr20\n\n";
    die;
} #end unless

open(LIST, "../$chr.phase_blocks.txt");

@blocks=();
@start=();
@end=();
@names=();

while (<LIST>) {
    chomp;
    @input=split(/\s+/,$_);
    push(@start,$input[1]);
    push(@end,$input[2]);
    push(@blocks,$input[3]);
} # end while

for ($i=0; $i<=$#start; $i++) {
    $phase_output = "$chr.$start[$i].$end[$i].out";
    $phase_input = "$chr.$start[$i].$end[$i].PHASE.gz";

    unless (-e $phase_output && -s $phase_output) {
	print STDERR "Error: File does not exist: $phase_output\n";
	die;
    } #end if

    unless (-e $phase_input && -s $phase_input) {
	print STDERR "Error: File does not exist: $phase_input\n";
	die;
    } #end if

    open(PHASE, $phase_output);
    $start_haps=0;
    $start_ind=0;
    %hap_set=();

    if ($i==0) {
	%haplotypes=();
	$output = "$chr.phased_block.$blocks[$i].out";
	open(OUTPUT,">$output");
	$block_count=1;
	$error_count=0;
	$new_block=1;
    } elsif ($blocks[$i]!=$blocks[$i-1]) { #last region in current block
	$new_block=1;
	$last_block=0;
	$num_ind = ($#names+1)/2;
	$num_sites=$#positions+1;
	$positions_as_scalar=join(" ",@positions);
	print OUTPUT "$num_ind\n$num_sites\nP $positions_as_scalar\n";
	for ($h=0; $h<=$#names; $h++) {
	    $ind=$names[$h];
	    print OUTPUT "#$ind\n$haplotypes{$ind}\n";
	    $length=length($haplotypes{$ind});
	} #end for
	$error_rate=$error_count/($block_count*$num_ind);
	$single_block=0;
	if ($num_sites!=$length) {
	    print STDERR "Block: $blocks[$i-1]; Length: $length; Num blocks: $block_count; Num sites: $num_sites; Error count: $error_rate\n";
	    print STDERR "Error:Number of positions in input does not match number of sites in output $phase_output!\n";
	    die;
	} #end if

	%haplotypes=();
	@positions=();
	$output = "$chr.phased_block.$blocks[$i].out";
	$block_count=1;
	$error_count=0;
	open(OUTPUT,">$output");
    } else {
	$block_count++;
    } #end if

    if ($blocks[$i+1]!=$blocks[$i]) { #next region starts new block
	$last_block=1;
    } #end if

    $positions=`zcat $phase_input | awk 'NR==3' | sed -e 's/P //'`;
    @dummy=split(/\s+/,$positions);    

    if ($new_block==1 && $last_block==1) { #only one region in block
	@positions=@dummy;
	$single_block=1;
#	print STDERR "Block: $blocks[$i-1]; Length: $length; Num blocks: $block_count; Num sites: $num_sites; Error count: $error_rate\n";
    } elsif ($i==0 || $new_block==1) { #remove last 50 sites 
	for ($p=0; $p<50;$p++) {
	    pop @dummy;
	} #end for
	push(@positions,@dummy);
    } elsif ($last_block==1) { #remove first 50 sites
	for ($p=0; $p<50;$p++) {
	    shift @dummy;
	} #end for
	push(@positions,@dummy);
    } else { #remove first and last 50 sites
	for ($p=0; $p<50;$p++) {
	    shift @dummy;
	    pop @dummy;
	} #end for
	push(@positions,@dummy);
    } #end else
#    if ($i<25 && $i>=1) {print STDERR "Positions: $#positions; $positions[0]-$positions[$#positions]\n";}

    while (<PHASE>) {
	chomp; 
	if ($_=~/BEGIN LIST_SUMMARY/) {
	    $start_haps=1;
	    next;
	} elsif ($_=~/END LIST_SUMMARY/) {
	    $start_haps=0;
	    next;
	} #end if
	if ($_=~/BEGIN BESTPAIRS_SUMMARY/) {
	    $start_ind=1;
	    next;
	} elsif ($_=~/END BESTPAIRS_SUMMARY/) {
	    $start_ind=0;
	    $new_block=0;
	    next;
	} #end elsif

	if ($start_haps==1) {
	    @hap_line=split(/\s+/,$_);
	    $hap_set{$hap_line[1]}=$hap_line[2];
#	    print STDERR "$hap_line[1]: $hap_line[2]\n";
	} elsif ($start_ind==1) {
	    @dum=split(/\s+/,$_);
	    $dum[0]=~s/(:|#)//g;
	    $ind=$dum[0];
	    $dum[1]=~s/(\(|\))//g;
	    @ind_haps=split(",",$dum[1]);
	    #extract 400 SNP sequence for current individual from haplotype hash
	    $ind_hap1=$hap_set{$ind_haps[0]};
	    $ind_hap2=$hap_set{$ind_haps[1]};
	    $hap1_name=$ind . "-1";
	    $hap2_name=$ind . "-2";

	    if ($i==0) {
		push(@names,$hap1_name,$hap2_name);
	    } #end if

	    if ($single_block==1) {
		$haplotypes{$hap1_name}=$ind_hap1;
		$haplotypes{$hap2_name}=$ind_hap2;
		next;
	    } elsif ($i==0 || $new_block==1) {
#		print STDERR "Ind: $ind; Hap 1: $ind_haps[0]; Hap 2: $ind_haps[1]\n";
		$haplotypes{$hap1_name}=substr($ind_hap1,0,350);
		$haplotypes{$hap2_name}=substr($ind_hap2,0,350);
		$last_haps{$hap1_name}=substr($ind_hap1,350,50); #351-400 from hap1,B
		$last_haps{$hap2_name}=substr($ind_hap2,350,50); #351-400 from hap2,D
		$second_last_haps{$hap1_name}=substr($ind_hap1,300,50); #301-350 from hap1,A
		$second_last_haps{$hap2_name}=substr($ind_hap2,300,50); #301-350 from hap2,C
		next;
#	    } elsif ($last_block==1) { #only one block in region
	    } #end if

	    $first_50_1=substr($ind_hap1,0,50); #1-50 from hap1,E
	    $first_50_2=substr($ind_hap2,0,50); #1-50 from hap2,G
	    $second_50_1=substr($ind_hap1,50,50); #51-100 from hap1,F
	    $second_50_2=substr($ind_hap2,50,50); #51-100 from hap2,H
	    
	    #calculate hamming distance for each haplotype combination, AB/CD and EF/GH
	    #A (sites 301-350 from B1) and B (sites 351-400 from B1), C (sites 301-350 from B1) and D (sites 351-400 from B1)
	    #E (sites 301-350 from B2) and F (sites 351-400 from B2), G (sites 301-350 from B2) and H (sites 351-400 from B2)
	    $hd_B_F=hd($last_haps{$hap1_name},$second_50_1);
	    $hd_C_G=hd($second_last_haps{$hap2_name},$first_50_2);
	    $hd_B_H=hd($last_haps{$hap1_name},$second_50_2);
	    $hd_C_E=hd($second_last_haps{$hap2_name},$first_50_1);
	    $x=$hd_B_F + $hd_C_G;  #Calculate the Hamming distance between B and F.  Add the Hamming distance between C and G.  Call this x.
	    $y=$hd_B_H + $hd_C_E;  #Calculate the Hamming distance between B and H.  Add the Hamming distance between C and E.  Call this y.

	    if ($last_block==1) { #keep all but first 50 SNPs for last block
		$hap1=substr($ind_hap1,50);
		$hap2=substr($ind_hap2,50);
		$test_length=length($hap1);
#		print STDERR "Block: $blocks[$i]; I: $i; Length=$test_length\n";
	    } else {   #only keep non-overlapping 300 SNPs, keeping 50 on either end for all other blocks
		$hap1=substr($ind_hap1,50,300);
		$hap2=substr($ind_hap2,50,300);
		$test_length=length($hap1);
		if ($test_length!=300) {
		    print STDERR "Error: region should only have 300 sites and instead has $test_length sites, output: $phase_output!\n";
		} #end if
	    } #end if

	    if ($x < $y) {  #If x<y, the the haplotypes are AF/CH.  Otherwise they are AH/CF.
		$haplotypes{$hap1_name}.=$hap1;
		$haplotypes{$hap2_name}.=$hap2;
#		if ($i==1) {print STDERR "Block: $blocks[$i]; I: $i; Ind: $ind; X: $x; Y: $y; BF: $hd_B_F; CG: $hd_C_G; BH: $hd_B_H; CE: $hd_C_E\n";}
	    } elsif ($x > $y) {  
#		if ($i==1) {print STDERR "Block: $blocks[$i]; I: $i; Ind: $ind; X: $x; Y: $y; BF: $hd_B_F; CG: $hd_C_G; BH: $hd_B_H; CE: $hd_C_E\n";}
		$haplotypes{$hap1_name}.=$hap2;
		$haplotypes{$hap2_name}.=$hap1;
	    } else {
#		print STDERR "Block: $blocks[$i]; I: $i; Ind: $ind; X: $x; Y: $y; BF: $hd_B_F; CG: $hd_C_G; BH: $hd_B_H; CE: $hd_C_E\n";
#		if ($i==1) {print STDERR "B:$last_haps{$hap1_name}\nD:$last_haps{$hap2_name}\nF:$second_50_1\nH:$second_50_2\n\nA:$second_last_haps{$hap1_name}\nC:$second_last_haps{$hap2_name}\nE:$first_50_1\nG:$first_50_2\n";}
#		print STDERR "Error: first 100 SNPs of current block and last 100 SNPs of previous block don't match\n";
		$haplotypes{$hap1_name}.=$hap1;
		$haplotypes{$hap2_name}.=$hap2;
		$error_count++;
	    } #end else

	    $last_haps{$hap1_name}=substr($ind_hap1,350,50); #351-400 from hap1,B
	    $last_haps{$hap2_name}=substr($ind_hap2,350,50); #351-400 from hap2,D
	    $second_last_haps{$hap1_name}=substr($ind_hap1,300,50); #301-350 from hap1,A
	    $second_last_haps{$hap2_name}=substr($ind_hap2,300,50); #301-350 from hap2,C
#	    $B_length=length( $last_haps{$hap1_name});
#	    $C_length=length($second_last_haps{$hap2_name});
#	    if ($B_length!=50 || $C_length!=50) {
#		print STDERR "Error: incorrect number of sites kept for comparison: B:$B_length, C:$C_length\n";
#	    } #end if
	} #end elsif
    } #end while
} #end for

$num_ind = ($#names+1)/2;
$num_sites=$#positions+1;
$positions_as_scalar=join(" ",@positions);
print OUTPUT "$num_ind\n$num_sites\nP $positions_as_scalar\n";

for ($h=0; $h<=$#names; $h++) {
    $ind=$names[$h];
    print OUTPUT "#$ind\n$haplotypes{$ind}\n";
    $length=length($haplotypes{$ind});
} #end for
$error_rate=$error_count/($block_count*$num_ind);
print STDERR "Block: $blocks[$i-1]; Length: $length; Num blocks: $block_count; Num sites: $num_sites; Error count: $error_rate\n";

#sub hd{ length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }

sub hd { #calculates hamming distance between two strings
    $a=$_[0];
    $b=$_[1];
    $ham_dist=0;
    $length_a=length($a);
    $length_b=length($b);

    if ($length_a!=$length_b) {
	warn "Error: strings of different lengths; $length_a, $length_b!";
    } #end if

    @array_a=split("", $a);
    @array_b=split("", $b);

    for ($t=0; $t<=$#array_a;$t++) {
	if ($array_a[$t] ne $array_b[$t]) {
	    $ham_dist++;
	} #end if
    } #end for

    return $ham_dist;

} #end hd

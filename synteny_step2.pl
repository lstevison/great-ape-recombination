#! /usr/bin/perl

#July 30, 2012
#Modified February 8, 2013
#program reads in output of synteny.pl
#concatenates consecutive synteny blocks with less than 15 SNPs skipped

$input = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Please provide input and output filenames on command line\n\n";
    die;
} #end unless

$output_two = "$output" . "2";

open(INPUT, $input);
open(OUTPUT, ">$output");
open(OUTPUT2, ">$output_two");

print OUTPUT2 "Region\t\#snps\tGorilla_chr\tGorilla_start\tGorilla_end\tGorilla_size\tHuman_chr\tHuman_start\tHuman_end\tHuman_size\tOrientation\n";
$header = (<INPUT>);

print STDERR "Reading output of synteny step 1...";

$region = 1;
$region_plus_one = 2;
@synteny_blocks = ();
@counts = ();
$count=0;
push(@counts, $count);

while (<INPUT>) {
    chomp;
    push(@synteny_blocks, $_);
    @input = split(/\t/, $_);
    if ($region == $input[0] || $input[0]=~/\-/) {
	$count++;
	next;
    } elsif ($input[0]==$region_plus_one) {
	push(@counts, $count);
	$count++;
	$region = $input[0];
	$region_plus_one = $region + 1;
    } else {
	print STDERR "$_\n";
    } #end if
} #end while

$count--;
push(@counts, $count);

print STDERR "done.\n";

$snp_count = 0;
$bad_snp_count = 0;
$region = 1;
@new_block = ();

for ($o=0; $o<$#counts; $o++) {
    $bad_snp_count=0;
    $current = $counts[$o];
    $next = $counts[$o+1];
    $second_next = $counts[$o+2];
    $third_next = $counts[$o+3];
    $start_current = $synteny_blocks[$current];
    @current_line = split(/\t/, $start_current);
    $start_next = $synteny_blocks[$next];
    @next_line = split(/\t/, $start_next);
    $start_second_next = $synteny_blocks[$second_next];
    @second_next_line = split(/\t/, $start_second_next);
    $size_of_middle = abs($second_next - $next + 1);
    $start_third_next = $synteny_blocks[$third_next];
    @third_next_line = split(/\t/, $start_third_next);
    $size_of_large_middle = abs($third_next - $next + 1);

#    print STDERR "Region: $current_line[0]-$next_line[0]; Start: $current_line[2], $next_line[2]; counts: $current-$next; Orientation: $current_line[5]/$next_line[5]; Size of middle: $size_of_middle\n";

    if ($snp_count==0) {
	$region_orientation=$current_line[5];
    } elsif ($snp_count>=1 && $region_orientation ne $current_line[5]) {
#print STDERR "Region: $current_line[0]-$next_line[0]; Start: $current_line[2], $next_line[2]; counts: $current-$next; Orientation: $current_line[5]/$next_line[5]; Size of middle: $size_of_middle\n";
#print STDERR "Orientation: $region_orientation; Size of block: $#new_block; SNP count: $snp_count; Current ori: $current_line[5]; Next_ori: $next_line[5]; 2nd next ori: $second_next_line[5]; 3rd next ori: $third_next_line[5]\n";
	if ($#new_block>=300) {
	    @first_line_block = split(/\t/, $new_block[0]);
	    $first_gor_pos = $first_line_block[2];
	    $first_hg18_pos = $first_line_block[4];
	    $block_gor_size = abs($last_gor_pos - $first_line_block[2]);
	    $block_hg18_size = abs($last_hg18_pos - $first_line_block[4]);
	    print OUTPUT2 "$region\t$snp_count\t$current_line[1]\t$first_line_block[2]\t$last_gor_pos\t$block_gor_size\t$current_line[3]\t$first_line_block[4]\t$last_hg18_pos\t$block_hg18_size\t$region_orientation\n";
	    for ($k=0; $k<=$#new_block; $k++) {
		print OUTPUT "$new_block[$k]\n";
	    }#end for
	    
	    $region++;
	    @new_block = ();
	    $bad_snp_count = 0;
	    $snp_count = 0;
	} else {
	    @new_block = ();
	    $bad_snp_count = 0;
	    $snp_count = 0;
	} #end elsif
	$region_orientation=$current_line[5];
    } #end elsif

    if ($current_line[5] eq $next_line[5] && $current_line[5] eq $region_orientation) { #orientation same between two blocks

#    print STDERR "Region: $current_line[0]-$next_line[0]; Start: $current_line[2], $next_line[2]; counts: $current-$next; Orientation: $current_line[5]/$next_line[5]; Size of middle: $size_of_middle\n";
	for ($i=$current; $i<$next; $i++) {
	    @input_array = split(/\t/, $synteny_blocks[$i]);
#	    if($input_array[0]==324 && $i==$current && $current_line[5]=~/plus/) {
#		$dist = $input_array[2]-$last_gor_pos;
#		print STDERR "Last gor pos: $last_gor_pos; current gor pos: $input_array[2]; Distance: $dist; SNP count: $snp_count; Last human pos: $last_hg18_pos; Current human pos: $input_array[4]; Ori: $current_line[5]\n";
#	    } #end if
	    if ($input_array[0]=~/\-/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/plus/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]<$next_line[4]) && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]<$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]<$last_gor_pos || $input_array[4]<$last_hg18_pos) && $snp_count>=1 && $current_line[5]=~/plus/) {
		$bad_snp_count++;
#		print STDERR "bad snp should find me: $input_array[2]!\n\n";
	    } elsif (($input_array[2]<$last_gor_pos || $input_array[4]>$last_hg18_pos) && $snp_count>=1 && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
		print STDERR "bad snp should find me: $input_array[2]!\n\n";
	    } elsif (($input_array[2]>$last_gor_pos || $input_array[4]<$last_hg18_pos) && $snp_count>=1 && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
		print STDERR "bad snp should find me: $input_array[2]!\n\n";
	    } else {
		$new_line = "$region\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[0]";
		push(@new_block, $new_line);
#		print OUTPUT "$new_line\n";
		$last_hg18_pos = $input_array[4];
		$last_gor_pos = $input_array[2];
		$region_orientation = $current_line[5];
		$snp_count++;
	    } #end else
	} #end for

	if ($next==$counts[$#counts]) {
	    @input_array = split(/\t/, $synteny_blocks[$next]);
	    $new_line = "$region\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[0]";
	    push(@new_block, $new_line);
	    $last_hg18_pos = $input_array[4];
	    $last_gor_pos = $input_array[2];
	    $snp_count++;
	} #end if

	$hg18_distance = abs($next_line[4] - $last_hg18_pos);
	$gor_distance = abs($next_line[2] - $last_gor_pos);

	if ($hg18_distance>=50000 || $gor_distance>=50000 || $next==$counts[$#counts]) {

	    if ($#new_block>=300) {
		@first_line_block = split(/\t/, $new_block[0]);
		$first_gor_pos = $first_line_block[2];
		$first_hg18_pos = $first_line_block[4];
		$block_gor_size = abs($last_gor_pos - $first_line_block[2]);
		$block_hg18_size = abs($last_hg18_pos - $first_line_block[4]);
		print OUTPUT2 "$region\t$snp_count\t$current_line[1]\t$first_line_block[2]\t$last_gor_pos\t$block_gor_size\t$current_line[3]\t$first_line_block[4]\t$last_hg18_pos\t$block_hg18_size\t$current_line[5]\n";
		for ($k=0; $k<=$#new_block; $k++) {
		    print OUTPUT "$new_block[$k]\n";
		}#end for

		$region++;
		@new_block = ();
		$bad_snp_count = 0;
		$snp_count = 0;
	    } else {
		@new_block = ();
		$bad_snp_count = 0;
		$snp_count = 0;
	    } #end elsif
	} #end if

    } elsif ($current_line[5] eq $second_next_line[5] && $region_orientation eq $current_line[5]) { #orientation same between current block and one after next and current region
	for ($i=$current; $i<$next; $i++) {
	    @input_array = split(/\t/, $synteny_blocks[$i]);
	    if ($input_array[0]=~/\-/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/plus/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]<$next_line[4]) && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]<$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]<$last_gor_pos || $input_array[4]<$last_hg18_pos) && $current_line[5]=~/plus/) {
		$bad_snp_count++;
#		print STDERR "bad snp should find me: $input_array[2]!\n\n";
	    } elsif ($snp_count>=1 && ($input_array[2]<$last_gor_pos || $input_array[4]>$last_hg18_pos) && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]>$last_gor_pos || $input_array[4]<$last_hg18_pos) && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
	    } else {
		$new_line = "$region\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[0]";
		push(@new_block, $new_line);
#		print OUTPUT "$new_line\n";
		$last_hg18_pos = $input_array[4];
		$last_gor_pos = $input_array[2];
		$region_orientation = $current_line[5];
		$snp_count++;
	    } #end else
	} #end for

	for ($i=$next; $i<$second_next; $i++) {
	    $bad_snp_count++;
	} #end for

	$hg18_distance = abs($second_next_line[4] - $last_hg18_pos);
	$gor_distance = abs($second_next_line[2] - $last_gor_pos);

	if ($hg18_distance>=50000 || $gor_distance>=50000 || $size_of_middle>300) {
	    if ($#new_block>=300) {
		@first_line_block = split(/\t/, $new_block[0]);
		$first_gor_pos = $first_line_block[2];
		$first_hg18_pos = $first_line_block[4];
		$block_gor_size = abs($last_gor_pos - $first_line_block[2]);
		$block_hg18_size = abs($last_hg18_pos - $first_line_block[4]);
		print OUTPUT2 "$region\t$snp_count\t$current_line[1]\t$first_line_block[2]\t$last_gor_pos\t$block_gor_size\t$current_line[3]\t$first_line_block[4]\t$last_hg18_pos\t$block_hg18_size\t$current_line[5]\n";
#		print STDERR "Region: $region; Length of block array: $#new_block; \#snps: $snp_count; Bad snp count: $bad_snp_count; Hg18 pos: $current_line[4]-$last_hg18_pos; GorGor pos: $current_line[2]-$last_gor_pos.\n";
		for ($k=0; $k<=$#new_block; $k++) {
		    print OUTPUT "$new_block[$k]\n";
		}#end for
		$region++;
		$bad_snp_count = 0;
		$snp_count = 0;
		@new_block = ();
	    } else {
		@new_block = ();
		$bad_snp_count = 0;
		$snp_count = 0;
	    } #end elsif
	} #end if

	$o++;

    } elsif ($current_line[5] eq $third_next_line[5] && $region_orientation eq $current_line[5]) { #orientation same between current block and one after next and current region
	for ($i=$current; $i<$next; $i++) {
	    @input_array = split(/\t/, $synteny_blocks[$i]);
	    if ($input_array[0]=~/\-/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/plus/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]<$next_line[4]) && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
	    } elsif (($input_array[2]<$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]<$last_gor_pos || $input_array[4]<$last_hg18_pos) && $current_line[5]=~/plus/) {
		print STDERR "bad snp should find me: $input_array[2]!\n\n";
		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]<$last_gor_pos || $input_array[4]>$last_hg18_pos) && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]>$last_gor_pos || $input_array[4]<$last_hg18_pos) && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
	    } else {
		$new_line = "$region\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[0]";
		push(@new_block, $new_line);
#		print OUTPUT "$new_line\n";
		$last_hg18_pos = $input_array[4];
		$last_gor_pos = $input_array[2];
		$region_orientation = $current_line[5];
		$snp_count++;
	    } #end else
	} #end for

	for ($i=$next; $i<$third_next; $i++) {
	    $bad_snp_count++;
	} #end for

	$hg18_distance = abs($third_next_line[4] - $last_hg18_pos);
	$gor_distance = abs($third_next_line[2] - $last_gor_pos);

	if ($hg18_distance>=50000 || $gor_distance>=50000 || $size_of_large_middle>300) {
	    if ($#new_block>=300) {
		@first_line_block = split(/\t/, $new_block[0]);
		$first_gor_pos = $first_line_block[2];
		$first_hg18_pos = $first_line_block[4];
		$block_gor_size = abs($last_gor_pos - $first_line_block[2]);
		$block_hg18_size = abs($last_hg18_pos - $first_line_block[4]);
		print OUTPUT2 "$region\t$snp_count\t$current_line[1]\t$first_line_block[2]\t$last_gor_pos\t$block_gor_size\t$current_line[3]\t$first_line_block[4]\t$last_hg18_pos\t$block_hg18_size\t$current_line[5]\n";
#		print STDERR "Region: $region; Length of block array: $#new_block; \#snps: $snp_count; Bad snp count: $bad_snp_count; Hg18 pos: $current_line[4]-$last_hg18_pos; GorGor pos: $current_line[2]-$last_gor_pos.\n";
		for ($k=0; $k<=$#new_block; $k++) {
		    print OUTPUT "$new_block[$k]\n";
		}#end for
		$region++;
		$bad_snp_count = 0;
		$snp_count = 0;
		@new_block = ();
	    } else {
		@new_block = ();
		$bad_snp_count = 0;
		$snp_count = 0;
	    } #end elsif
	} #end if

	$o++;
	$o++;
    } else {
	for ($i=$current; $i<$next; $i++) {
	    @input_array = split(/\t/, $synteny_blocks[$i]);
	    if ($input_array[0]=~/\-/) {
		$bad_snp_count++;
#	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/plus/) {
#		$bad_snp_count++;
#	    } elsif (($input_array[2]>$next_line[2] || $input_array[4]<$next_line[4]) && $current_line[5]=~/minus1/) {
#		$bad_snp_count++;
#	    } elsif (($input_array[2]<$next_line[2] || $input_array[4]>$next_line[4]) && $current_line[5]=~/minus2/) {
#		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]<$last_gor_pos || $input_array[4]<$last_hg18_pos) && $current_line[5]=~/plus/) {
		$bad_snp_count++;
		print STDERR "bad snp should find me: $input_array[2]!\n\n";
	    } elsif ($snp_count>=1 && ($input_array[2]<$last_gor_pos || $input_array[4]>$last_hg18_pos) && $current_line[5]=~/minus1/) {
		$bad_snp_count++;
	    } elsif ($snp_count>=1 && ($input_array[2]>$last_gor_pos || $input_array[4]<$last_hg18_pos) && $current_line[5]=~/minus2/) {
		$bad_snp_count++;
	    } else {
		$new_line = "$region\t$input_array[1]\t$input_array[2]\t$input_array[3]\t$input_array[4]\t$input_array[5]\t$input_array[0]";
		push(@new_block, $new_line);
#		print OUTPUT "$new_line\n";
		$last_hg18_pos = $input_array[4];
		$last_gor_pos = $input_array[2];
		$region_orientation = $current_line[5];
		$snp_count++;
	    } #end else
	} #end for
    } #end elsif

} #end for


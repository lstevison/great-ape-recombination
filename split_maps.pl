#! /usr/bin/perl

#LSS October 22, 2012
#Modified October 3, 2014

#program reads in lifted_maps_min_chunk_5000_SNPs_50kb_gaps_10CEU_10YRI_with_chimp_super_clean_genotype_map_and_HapMap_Pops.txt and splits it into corresponding files for each map and converts rec rate to rho/kb for all

$input = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_10CEU_10YRI_with_chimp_super_clean_genotype_map_and_HapMap_Pops.txt";

open(INPUT, $input);
$header = (<INPUT>);
$panmap = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_chimp_genotype_map_panTro2_coords.txt";
#$panmap1 = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_chimp_super_clean_genotype_map.txt";
#$panmap2 = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_chimp.txt";
#$CEU_map = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_CEU_HapMap_Pops.txt";
#$YRI_map = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_YRI_HapMap_Pops.txt";
#$hapmap = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_HapMap_Pops.txt";
#$TEN_CEU = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_10CEU_Pops.txt";
#$TEN_YRI = "lifted_maps_min_chunk_5000_SNPs_50kb_gaps_10YRI_Pops.txt";

open(PANMAP, ">$panmap");
#open(PANMAP1, ">$panmap1");
#open(PANMAP2, ">$panmap2");
#open(CEU, ">$CEU_map");
#open(YRI, ">$YRI_map");
#open(HAPMAP, ">$hapmap");
#open(TEN_YRI, ">$TEN_YRI");
#open(TEN_CEU, ">$TEN_CEU");

#open(BLOCKS, ">Panmap.syntenic_blocks.hg18.BED");

$first_interval = 1;

while (<INPUT>) {
    chomp; 
    @input_array = split(/\t/, $_);
    $chr = $input_array[5];
    $region = $input_array[0];
    $start = $input_array[6];
    $mb_start = $start/1000000;

    $pantro_chr = $input_array[1];
    $pantro_start = $input_array[2];
    $mb_pantro_start = $pantro_start/1000000;


    if ($first_interval==1) {
	$last_region = $region;
	$last_start = $input_array[6];
	$last_mb_start = $mb_start;
	$last_chimp_rate2 = $input_array[3];
	$last_chimp_rate1 = $input_array[4];
	$last_hapmap_rate = $input_array[7];
	$last_CEU_rate = $input_array[10];
	$last_YRI_rate = $input_array[11];
	$last_TEN_CEU_rate = $input_array[8];
	$last_TEN_YRI_rate = $input_array[9];
	$last_pantro_start = $input_array[2];
	$last_mb_pantro_start = $mb_pantro_start;
#	print BLOCKS "Interval\tChr\tStart\tEnd\n";
#	print BLOCKS "$input_array[0]\t$input_array[5]\t$input_array[6]\t";
	print STDERR "Now processing interval $input_array[0]...";
	$first_interval = 0;
	next;
    } #end if

    if ($region != $last_region) {
#	print BLOCKS "$last_start\n$input_array[0]\t$input_array[5]\t$input_array[6]\t";
	print STDERR "$input_array[0]...";
	$last_region = $region;
	$last_start = $input_array[6];
	$last_mb_start = $mb_start;
	$last_chimp_rate2 = $input_array[3];
	$last_chimp_rate1 = $input_array[4];
	$last_hapmap_rate = $input_array[7];
	$last_CEU_rate = $input_array[10];
	$last_YRI_rate = $input_array[11];
	$last_TEN_CEU_rate = $input_array[8];
	$last_TEN_YRI_rate = $input_array[9];
	$last_pantro_start = $input_array[2];
	$last_mb_pantro_start = $mb_pantro_start;
	next;
    } elsif ($region==$last_region) {

	$kb_distance = (abs($start - $last_start))/1000; 
	$pantro_kb_distance = (abs($pantro_start - $last_pantro_start))/1000; 

	$chimp_rate = (abs($input_array[3] - $last_chimp_rate2))/$pantro_kb_distance;

	$chimp_rate1 = (abs($input_array[4] - $last_chimp_rate1))/$kb_distance;
	$chimp_rate2 = (abs($input_array[3] - $last_chimp_rate2))/$kb_distance;
	$hapmap_rate = (((abs($input_array[7] - $last_hapmap_rate))/100)*(4*10040))/$kb_distance;
	$CEU_rate = (((abs($input_array[10] - $last_CEU_rate))/100)*(4*10040))/$kb_distance;
	$YRI_rate = (((abs($input_array[11] - $last_YRI_rate))/100)*(4*19064))/$kb_distance;
	$TEN_CEU_rate = (abs($input_array[8] - $last_TEN_CEU_rate))/$kb_distance;
	$TEN_YRI_rate = (abs($input_array[9] - $last_TEN_YRI_rate))/$kb_distance;
	
	print PANMAP "$pantro_chr\t$last_mb_pantro_start\t$mb_pantro_start\t$chimp_rate\t$region\n";
#	print PANMAP1 "$chr\t$last_mb_start\t$mb_start\t$chimp_rate1\t$region\n";
#	print PANMAP2 "$chr\t$last_mb_start\t$mb_start\t$chimp_rate2\t$region\n";
#	print HAPMAP "$chr\t$last_mb_start\t$mb_start\t$hapmap_rate\n";
#	print CEU "$chr\t$last_mb_start\t$mb_start\t$CEU_rate\n";
#	print YRI "$chr\t$last_mb_start\t$mb_start\t$YRI_rate\n";
#	print TEN_CEU "$chr\t$last_mb_start\t$mb_start\t$TEN_CEU_rate\n";
#	print TEN_YRI "$chr\t$last_mb_start\t$mb_start\t$TEN_YRI_rate\n";
	
	$last_region = $region;
	$last_start = $input_array[6];
	$last_mb_start = $mb_start;
	$last_chimp_rate2 = $input_array[3];
	$last_chimp_rate1 = $input_array[4];
	$last_hapmap_rate = $input_array[7];
	$last_CEU_rate = $input_array[10];
	$last_YRI_rate = $input_array[11];
	$last_TEN_CEU_rate = $input_array[8];
	$last_TEN_YRI_rate = $input_array[9];
	$last_pantro_start = $input_array[2];
	$last_mb_pantro_start = $mb_pantro_start;
    } #end elsif

} #end while

#print BLOCKS "$input_array[6]\n";

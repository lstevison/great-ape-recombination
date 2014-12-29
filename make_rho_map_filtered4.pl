#! /usr/bin/perl

#Last modified February 11, 2013
#Created: November 28, 2011                                                                                                      

#for a series of LDhat outputs, this program combines the locs file and res file to map recombination rate
# it also reads in the consum file to plot a summary of MCMC convergence for each run

$chromosome = $ARGV[0];
$path_to_inputs = "./LDhat-inputs/$chromosome";
$path_to_outputs = "./LDhat-outputs/$chromosome";

open(LDHAT_LIST, "$path_to_outputs/$chromosome.input-list_new.txt");
open(MASTER_OUTPUT, ">$path_to_outputs/$chromosome-map.txt");

print STDERR "Start coord\tend coord\tNumber SNPs\tNumber SNPs high rho(filter1)\tNumber SNPs set to zero(filter2)\n";

while (<LDHAT_LIST>) {
    chomp;
    $ldhat_input = $_; 
    @header = split(/\./, $ldhat_input);
    chomp($header[2]);
    $size_of_region = ($header[2] - $header[1]) + 1;
    push(@ldhat_file, $ldhat_input);
    push(@region_sizes, $size_of_region);
} #end while

open(R_SCRIPT, ">maps-$chromosome.r");    					#open R-script

for ($z=0; $z<=$#ldhat_file; $z++) {

    @rho = ();
    @snp_loc_mb = ();
    @snp_loc_kb = ();
    @U95 = ();
    @L95 = ();

    open(LOCS_FILE, "$path_to_inputs/$ldhat_file[$z].ldhat.locs");     		#gives bp coordinate for each marker, x-axis
    open(RES_FILE, "$path_to_outputs/res.$ldhat_file[$z].txt");	       		#gives rho estimate for each interval, y-axis
    $first_time_in_loop = 1;
    $while_counter = 0;
    
    while (<LOCS_FILE>) {			      						#read in locs file
	chomp;
	$locs_input = $_; 
	$while_counter++;
	
	if ($first_time_in_loop == 1) {		       					#first line of locs file
	    @first_line = split(/\s+/, $locs_input);
	    $num_snps = $first_line[0];						 #capture number of snps for this file
	    push(@number_snps, $num_snps);
	    $max_coord = $first_line[1]/1000;
	    $max_kb_coord = $first_line[1];
	    $first_time_in_loop = 0;
	} elsif ($while_counter <= 100) {
	    next;
	} #end elsif

	if ($while_counter==101) {
	    $min_kb_coord = $locs_input;
	} #end if

    } #end while

    $min_rho=0;
    $max_rho=0;	
    $loop_counter=0;
    
    while (<RES_FILE>) {				 				#read in res file
	chomp;
	$res_input = $_; 
	$loop_counter++;
	
	if ($loop_counter > 100) {							#discard first two lines of res file
	    @line_item = split(/\t/, $res_input);
	    $snp_in_mb = $line_item[0]/1000;
	    push(@snp_loc_mb, $snp_in_mb);
	    $res_input2 = $line_item[0]/1;
	    push(@snp_loc_kb, $res_input2);

	    $mean_rho = $line_item[1]/1;		     					#capture mean 
	    push(@rho, $mean_rho);						        #capture rho value into an array
	    push(@L95,$line_item[3]);
	    push(@U95, $line_item[4]);
#	    $snp_in_mb = $line_item[0]/1000;				     	   #convert coordinate from kb to mb
#	    push(@snp_loc_mb, $snp_in_mb);				       	   #capture snp location into an array
#	    push(@snp_loc_kb, $line_item[0]);					   #capture snp location into an array
	} #end elsif
    } #end while

    if ($z==0) {
	$min_coord = $snp_loc_mb[0];
    } #end if

    @filter1_array=();

    for ($filter1=0; $filter1<=$#rho; $filter1++) {
	$interval_size = ($snp_loc_kb[$filter1+1]-$snp_loc_kb[$filter1]) + 0.001;
	$rho_per_kb = $rho[$filter1];	
	$actual_rho = $rho_per_kb*$interval_size;

	if ($actual_rho >= 100) {
	    push(@filter1_array,$snp_loc_kb[$filter1]); 
	} #end if

    } #end filter 1 loop

    @filter1_array = sort { $a <=> $b } (@filter1_array);	
    %filter1_hash = map { $_ => 1 } @filter1_array;
    @unique_filter1 = keys %filter1_hash;
    @unique_filter1 = sort { $a <=> $b } (@unique_filter1);	

    $length_filter1_array = $#unique_filter1 + 1;
    $j=0;
    @filter2_array=();

    if ($length_filter1_array>0) {
	
	for ($filter2=0; $filter2<=$#snp_loc_kb; $filter2++) {

	    if ($snp_loc_kb[$filter2]==$unique_filter1[$j] && $filter2<=50) {
		$j++;
		$start_for_loop = 0;
		$up25_snp_pos = $snp_loc_kb[$start_for_loop];
		$stop_for_loop = $filter2 + 50;
		$down25_snp_pos = $snp_loc_kb[$stop_for_loop];
	    } elsif ($snp_loc_kb[$filter2]==$unique_filter1[$j] && $filter2>50) {
		$j++;
		$start_for_loop = $filter2 - 50;
		$up25_snp_pos = $snp_loc_kb[$start_for_loop];
		$stop_for_loop = $filter2 + 50;
		$down25_snp_pos = $snp_loc_kb[$stop_for_loop];
	    } elsif ($j==$length_filter1_array) {
		last;
	    } else {
		next;
	    } #end if
	    
	    for ($n=$start_for_loop; $n<=$stop_for_loop; $n++) {
		if (defined($snp_loc_kb[$n])) {
		    push(@filter2_array,$snp_loc_kb[$n]);		
		} #end if
	    } #end n-for

	} #end filter 2 loop	

	@filter2_array = sort { $a <=> $b } (@filter2_array);	
	%filter2_hash = map { $_ => 1 } @filter2_array;
	#print STDERR "Filter 2 array made into hash!\n";
	@unique_filter2 = keys %filter2_hash;
	@unique_filter2 = sort { $a <=> $b } (@unique_filter2);	
	$no_filtered_sites = 0;

    } else {
	@unique_filter2 = ();
	$no_filtered_sites = 1;
    } #end else
	
    $m=0;
    $n=0;
    $num_unique_snps_filter2 = $#unique_filter2 + 1; 
    $stop_print_temp_file = $#snp_loc_kb - 100;
    $zero_counter = 0;
    print STDERR "$snp_loc_mb[0]\t$snp_loc_mb[$stop_print_temp_file]\t";
    
    for ($i=0; $i<$stop_print_temp_file; $i++) {
	$interval_midpoint = ($snp_loc_mb[$i+1] + $snp_loc_mb[$i])/2;
	$rho_per_kb = $rho[$i];
	$Upper_95 = $U95[$i];
	$Lower_95 = $L95[$i];

	if ($snp_loc_kb[$i]==$unique_filter2[$m] && $no_filtered_sites==0) {
	    if ($unique_filter1[$n]==$unique_filter2[$m]) {
		print MASTER_OUTPUT "$chromosome\t$snp_loc_mb[$i]\t$snp_loc_mb[$i+1]\t$interval_midpoint\t100\t100\t100\n";			#print into temp-locs file
		$zero_counter++;
		$m++;
		$n++;
		next;
	    } else {
		print MASTER_OUTPUT "$chromosome\t$snp_loc_mb[$i]\t$snp_loc_mb[$i+1]\t$interval_midpoint\t150\t150\t150\n";			#print into temp-locs file
		$zero_counter++;
		$m++;
		next;
	    } #end else
	} elsif ($m>$num_unique_snps_filter2 || $no_filtered_sites==1) {

	    if ($rho_per_kb > $ max_rho) {
		$max_rho = $rho_per_kb;					#re-assign max rho if current rho is greater
	    } elsif ($rho_per_kb < $min_rho) {
		$min_rho = $rho_per_kb;	
	    } #end elsif

	    print MASTER_OUTPUT "$chromosome\t$snp_loc_mb[$i]\t$snp_loc_mb[$i+1]\t$interval_midpoint\t$rho_per_kb\t$Lower_95\t$Upper_95\n";	     #print into temp-locs file
	    next;
	} else {

	    if ($rho_per_kb > $ max_rho) {
		$max_rho = $rho_per_kb;					#re-assign max rho if current rho is greater
	    } elsif ($rho_per_kb < $min_rho) {
		$min_rho = $rho_per_kb;	
	    } #end elsif

	    print MASTER_OUTPUT "$chromosome\t$snp_loc_mb[$i]\t$snp_loc_mb[$i+1]\t$interval_midpoint\t$rho_per_kb\t$Lower_95\t$Upper_95\n";		#print into temp-locs file
	    next;
	} #end else
	
    } #end i-for
    
    $num_snps_filter1 = $#unique_filter1+1;
    $num_snps_filter2 = $#unique_filter2+1;
    
    print STDERR "$stop_print_temp_file\t$num_snps_filter1\t$zero_counter\n";
    #for each output, read in consum file
    print R_SCRIPT "consum$z<-read.table(\"$path_to_outputs/consum.$ldhat_file[$z].txt\", header=FALSE, sep=\"\\t\")\n";
    
} #end z-for

#read in master output
print R_SCRIPT "main_map<-read.table(\"$path_to_outputs/$chromosome-map.txt\", header=FALSE, sep=\"\\t\")\n";

#make subset to remove high rho sites from plotting
print R_SCRIPT "new_map=subset(main_map, main_map\$V5 < 100)\n";
print R_SCRIPT "coord=new_map\$V4\n";
print R_SCRIPT "rate=new_map\$V5\n";

#open pdf to draw map
print R_SCRIPT "png(file=\"rho_figure_$chromosome.png\", height=700, width=1500, units=\"px\", pointsize=\"20\")\n";
print R_SCRIPT "par(mar=c(3,4,3,4)+0.1)\n";

#print blank map first and text for both axes
print R_SCRIPT "plot(c($min_coord,$max_coord), c(0, 5), type=\"n\", xlab=\"\", ylab=\"\", main=\"Broad-scale and fine-scale recombination for $chromosome\")\n";
print R_SCRIPT "mtext(\"Position (Mb)\", side=1, line=2, font=2)\n";
print R_SCRIPT "mtext(\"Fine-scale recombination rate - rho/kb\", side=4, line=3, cex=1.1)\n";
print R_SCRIPT "mtext(\"Broad-scale recombination rate - rho/kb\", side=2, line=3, cex=1.1)\n";

#right axis is fine-scale rate
print R_SCRIPT "par(new=TRUE)\n";
print R_SCRIPT "plot(c($min_coord,$max_coord), c(0, 50), type=\"n\", xaxt=\"n\",yaxt=\"n\", xlab=\"\", ylab=\"\", main=\"\")\n";
print R_SCRIPT "axis(4)\n";
print R_SCRIPT "lines(coord, rate, type=\"s\", lwd=1, col=\"grey\")\n";	      #plots full map

#left axis is broad scale, perform loess smoothing of map and plot
print R_SCRIPT "par(new=TRUE)\n";
print R_SCRIPT "plot(c($min_coord,$max_coord), c(0, 5), type=\"n\", xaxt=\"n\",yaxt=\"n\", xlab=\"\", ylab=\"\", main=\"\")\n";
print R_SCRIPT "axis(2)\n";
print R_SCRIPT "map.loess<-loess(rate~coord,span=0.02,data.frame(x=coord,y=rate))\n";
print R_SCRIPT "map.predict<-predict(map.loess,data.frame(x=coord))\n";
print R_SCRIPT "lines(coord, map.predict, type=\"s\", lwd=3, col=rgb\(0,0,0.5\))\n";	      #plots loess smoothed map

print R_SCRIPT "dev.off()\n";	 	#writes png to file

#plot convergence summary
print R_SCRIPT "pdf(file=\"consum_figure_$chromosome.pdf\", height=8.5, width=11)\n";		#open pdf to draw graphs
print R_SCRIPT "par(mfrow=c(2,2))\n";				#split for multiple graphs

for ($b=0; $b<=$#ldhat_file; $b++) {		#iteratively plots both convergence summary graphs
    print R_SCRIPT "plot(consum$b\$V2, ylab=\"Number of blocks\", xlab=\"\#MCMC iterations (x10000)\", type=\"l\")\n";
    print R_SCRIPT "title(main=\"Convergence Summary: $ldhat_file[$b]\")\n";
    print R_SCRIPT "plot(consum$b\$V1, ylab=\"Total map length\", xlab=\"\#MCMC iterations (x10000)\", type=\"l\")\n";
    print R_SCRIPT "title(main=\"Number of SNPs: $number_snps[$b]; Size of region: $region_sizes[$b] bp\")\n";
} #end b-for

print R_SCRIPT "dev.off()\n";			  #writes to pdf file

#print STDERR "Number regions with high rho = $high_rho_counter; Number SNPs with high rho: $num_snps_highrho\n";

system("R --vanilla <maps-$chromosome.r");			  #execute R script generated within program

#make hotspot usage plot
#make chromosome comparison plot

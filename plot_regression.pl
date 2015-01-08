#! /usr/bin/perl

#October 22, 2012
#Last modified: August 11, 2014

$input = $ARGV[0];
$output = $ARGV[1];

unless ($#ARGV==1) {
	print STDERR "Please specify on command line input and output filenames.\n";
	die;
} #end unless

open(INPUT, $input);
open(OUTPUT, ">$output");
open(PLOT, ">plot_regression.r");

@start_pos = ();
@end_pos = ();
@rec_rate = ();
$total_map_length = 0;
$total_physical_length = 0;
%hash = ();
$previous_recrate = 0;
$previous_dist = 0;
$header = (<INPUT>);

while (<INPUT>) {
	chomp;
	@input_array = split(/\s+/, $_);
	push(@start_pos,$input_array[1]);
	push(@end_pos,$input_array[2]);
	push(@rec_rate,$input_array[3]);
} #end while

for ($i=0; $i<$#start_pos; $i++) {
	$distance = abs($end_pos[$i] - $start_pos[$i]) + 1;
#	$kb_distance = $distance * 1000;
	#print STDERR "$distance\n";
#	$rho = $rec_rate[$i]*$kb_distance;
	$total_map_length += $rec_rate[$i];
	$total_physical_length += $distance;

	$key = "$distance" . ':' . "$start_pos[$i]" . ':' . "$rec_rate[$i]";
	$value = $rec_rate[$i];
	$hash{$key} = $value;
	#print STDERR "Value: $value; Key: $key\n";
} #end for

print STDERR "Total map length: $total_map_length; Total physical length: $total_physical_length\n";

#sort hash based on rec rate (ascending)

foreach $key (sort { $hash{$a} <=> $hash{$b} } keys %hash) {
	@base_array = split (/\:/, $key);
	#print STDERR "Rho: $base_array[2]; Distance: $base_array[0]\n";
	$prop_gen = (($base_array[2] / $total_map_length)*100) + $previous_recrate;
	$prop_dist = (($base_array[0] / $total_physical_length)*100) + $previous_dist;
	print OUTPUT "$prop_dist\t$prop_gen\n";
	$previous_recrate = $prop_gen;
	$previous_dist = $prop_dist;
}

#plot OUTPUT
print PLOT "CEU_human<-read.table(\"CEU_human.txt\", header=FALSE)\n";
print PLOT "YRI_human<-read.table(\"YRI_human.txt\", header=FALSE)\n";
print PLOT "clean.chimp<-read.table(\"chimp_clean.txt\", header=FALSE)\n";
print PLOT "png(file=\"prop_seq.png\", height=7.5, width=7, units=\"in\", pointsize=\"25\", bg= \"transparent\", res=320)\n";
print PLOT "par(mar=c(3,3,2,0.5))\n";
print PLOT "plot(c(0,100), c(0,100), type=\"n\", xlab=\"\", ylab=\"\", main=\"\")\n";
print PLOT "lines(YRI_human\$V1, YRI_human\$V2, type=\"s\", col=\"sienna\", lwd=7)\n";
print PLOT "mtext(\"Proportion of Recombination\", side=1, line=2)\n";
print PLOT "mtext(\"Proportion of Sequence\", side=2, line=2)\n";
print PLOT "mtext(\"Hotspot Usage\", side=3, line=1)\n";
print PLOT "lines(clean.chimp\$V1, clean.chimp\$V2, type=\"s\", col=\"green\", lwd=5)\n";
print PLOT "lines(CEU_human\$V1, CEU_human\$V2, type=\"s\", col=\"blue\", lwd=7)\n";
print PLOT "legend(\"topleft\", c(\"CEU HapMap\",\"Chimp\",\"YRI HapMap\"), col=c(\"blue\",\"green\",\"sienna\"), lwd=3, bty=\"n\")\n";

print PLOT "dev.off()\n";

#system ("R --vanilla <plot_regression.r");

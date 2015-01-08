#! /usr/bin/perl
 
#7/5/2008 LSS 
#Modified September 10, 2014
 
sub usage { 
    print STDERR "usage: extract sequences <Source sequence file> <list of positions> <Output filename>\n\n"; 
    print STDERR "Extracts positions from source sequence file and outputs a fasta delimited file.\n"; 
    print STDERR "Creator: Laurie Stevison\n"; 
    exit; 
} #end sub 
 
if ($#ARGV  != 2) { # zero means one argument 
    &usage; 
} #end if 

$source = $ARGV[0]; 
$positions = $ARGV[1]; 
$output = $ARGV[2]; 
$output2 = $output;
$output2=~s/FASTA/GCcount.txt/;

open(READSRC, $source); 
open(READPOS, $positions); 
open(WRITEOUT, ">$output"); 
open(OUT2, ">$output2"); 

print STDERR "Reading input file...";

@chr=();
@start=();
@end=();

while (<READPOS>) {
    chomp;
    @input_array = split(/\s+/, $_);
    push(@chr , $input_array[0]);
    push(@start , $input_array[1]);
    push(@end , $input_array[2]);
} #end while

print STDERR "done.\n\n";

print STDERR "Reading source file...";

while (<READSRC>) {
    chomp;
    if (/>/) { 
        $seq_name = $'; #'
            push(@grp_names, $seq_name);
        if ($newseq==1) {push(@sequences, $sequence);} 
        $sequence = "";
        $newseq=0;
        next; 
    } else {
        $sequence .= $_;
        $newseq=1; 
    } #end else
} #end while
push(@sequences, $sequence);

print "done.\n\n";

print STDERR "Extracting positions...\n";

for ($i = 0; $i <= $#start; $i++) {
    
    for ($k = 0; $k <= $#grp_names; $k++) {
        if ($grp_names[$k] eq $chr[$i]) {
            $hold_sequence = $sequences[$k];
        } else {
            next;
        } #end else
    } #end for
    
    $interval_size = $end[$i]-$start[$i]+1;
    $extract = substr($hold_sequence,  $start[$i], $interval_size);
    $a=($extract=~tr/Aa//);
    $t=($extract=~tr/Tt//);
    $g=($extract=~tr/Gg//);
    $c=($extract=~tr/Cc//);
    $total=$a+$t+$g+$c;
    if ($total != 0) {
	$gc= sprintf("%.2f", (($g+$c)/$total)*100 );
    } else {
	$gc=0;
    } #end else                                                                                                                                                                  

    $percent_N=sprintf("%.2f",(($interval_size-$total)/$interval_size)*100);

    print WRITEOUT ">$chr[$i].$start[$i]:$end[$i]\n";
    print WRITEOUT "$extract\n";
    print OUT2 "$chr[$i]\t$start[$i]\t$end[$i]\t$gc\t$percent_N\n";    
} #end for


print STDERR "done.\n\n";

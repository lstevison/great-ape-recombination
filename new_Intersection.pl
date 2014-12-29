#! /usr/bin/perl

#July 10, 2012

#program reads in output of awk intersection in kept.sites format and output of INFO field to convert all coordinates in intersection file back to NHP coordinates

$intersect = $ARGV[0];
$INFO = $ARGV[1];
$OUTPUT = $ARGV[2];

unless ($#ARGV==2) {
    print STDERR "Usage: Program takes in intersection output in hg18 coordinates and original INFO file with both NHP and hg18 coordinates to converts back to NHP coordinates. Please specify on command line the input filenames and an output filename.\n\n";
    die;
} #end unless

print STDERR "Running intersection to NHP coordinate program...";
open(HG18, $intersect);
open(INFO, $INFO);
open(OUTPUT, ">$OUTPUT");

#remove original headers
$header1 = (<HG18>);
$header2 = (<INFO>);
#replace with new header
print OUTPUT "CHROM\tPOS\n";
#@hg18 = ();
%info = ();
#@info_array = ();

while (<INFO>) {
    chomp;
    $info_line = $_;
    @site = split(/\t/, $info_line);
    $info{$site[4]} = "$site[0]\t$site[1]";
#    push (@info_array, $site[4]) ;
} #end while


while (<HG18>) {
    chomp;
    $line = $_;
    @whole_line = split(/\t/, $line);
    $HG18 = "$whole_line[0]:$whole_line[1]";
    if (defined($info{$HG18})) {
	print OUTPUT "$info{$HG18}\n";
    } #end if
#    push(@hg18, $HG18);
} #end while

print STDERR "done.\n";

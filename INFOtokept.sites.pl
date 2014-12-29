#! /usr/bin/perl

#July 2, 2012

#program reads in OUTPUT from INFO field extraction and converts to same format as kept.sites file (removing failed hg18 RL sites)

$INFO = $ARGV[0];
$OUTPUT = $ARGV[1];

unless ($#ARGV==1) {
    print STDERR "Usage: Program converts INFO output to kept.sites format. Please specify on command line the input filename.\n\n";
    die;
} #end unless

print STDERR "Running conversion program on INFO output...";
open(INFO, $INFO);
open(OUTPUT, ">$OUTPUT");

#remove original header
$header = (<INFO>);
#replace with new header
print OUTPUT "CHROM\tPOS\n";

while (<INFO>) {
    chomp;
    $line = $_;
    @site = split(/\t/, $line);
    if ($site[4]!~/failed/) {
	@hg18 = split(":", $site[4]);
	print OUTPUT "$hg18[0]\t$hg18[1]\n";
    } #end if
} #end while

print STDERR "\tdone.\n";

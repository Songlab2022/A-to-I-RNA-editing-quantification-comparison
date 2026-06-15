#!/usr/bin/env perl
############################################################
# Query editing level of sites in a BAM file
# Usage: perl Query_Editing_Level.pl <site_list.bed> <indexed.bam> <output.txt> [min_base_qual]
############################################################
use warnings;
use strict;
require "parse_pileup.pl";   

if (@ARGV < 3) {
    die "Usage: perl Query_Editing_Level.pl <site_bed> <bam_file> <output.txt> [min_base_qual]\n";
}
my ($inputfile, $bamfile, $outputfile, $minbasequal) = @ARGV;
$minbasequal = 30 unless defined $minbasequal;  #ngs=30；lrs=7

my $minmapqual = 20;
my $sampath = "samtools";         
my $genomepath = "";             
my $offset = 33;                  # Phred+33

die "Error: please set \$genomepath in the script\n" unless $genomepath;

my $bedtemp = $outputfile . '.bed';
system("awk \'{print \$1\"\t\"\$2-1\"\t\"\$2}\' $inputfile > $bedtemp");
my $piletemp = $outputfile . '.pileup';
system("$sampath mpileup -A -B -d 1000000 -f $genomepath -l $bedtemp $bamfile > $piletemp");

my %sitehash;
open my $PILEUP, "<", $piletemp or die "Cannot open $piletemp: $!";
while (<$PILEUP>) {
    chomp;
    my ($chr, $position, $refnuc, $coverage, $pile, $qual) = split;
    my $location = join '_', $chr, $position;
    my ($refnuccount, $acount, $tcount, $ccount, $gcount) = parse_pileup($_, $minbasequal, $offset);
    my $counts = join ',', $refnuccount, $ccount, $gcount;
    $sitehash{$location} = $counts;
}
close $PILEUP;
unlink $bedtemp, $piletemp;

open my $INPUT, "<", $inputfile or die "Cannot open $inputfile: $!";
open my $OUTPUT, ">", $outputfile or die "Cannot write $outputfile: $!";
print $OUTPUT "#chrom\tposition\tstrand\tcoverage\teditedreads\teditlevel\n";

while (<$INPUT>) {
    chomp;
    my @fields = split;
    my ($chr, $position, $strand) = ($fields[0], $fields[1], $fields[2]);
    my $location = join '_', $chr, $position;
    if ($sitehash{$location}) {
        my ($refcount, $ccount, $gcount) = split /,/, $sitehash{$location};
        my $newmismatch = ($strand eq '+') ? $gcount : $ccount;
        my $newcov = $refcount + $newmismatch;
        if ($newcov) {
            my $varfreq = sprintf("%.3f", $newmismatch / $newcov);
            print $OUTPUT "$fields[0]\t$fields[1]\t$strand\t$newcov\t$newmismatch\t$varfreq\n";
        } else {
            print $OUTPUT "$fields[0]\t$fields[1]\t$strand\t0\t0\tN/A\n";
        }
    } else {
        print $OUTPUT "$fields[0]\t$fields[1]\t$strand\t0\t0\tN/A\n";
    }
}
close $INPUT;
close $OUTPUT;

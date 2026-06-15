#!/usr/bin/env perl
use strict;
use Getopt::Long;
use File::Basename;

my ($in, $out, $ref, $sites, $help);
GetOptions(
    "in=s"   => \$in,
    "out=s"  => \$out,
    "ref=s"  => \$ref,
    "sites=s"=> \$sites,
    "help"   => \$help
);

if ($help) {
    print <<"EOF";
Usage: perl editing_site_call.pl --in <input_dir> --out <output_dir> --ref <genome.fa> --sites <sites.bed>

Description:
    Scan all .bam files in the input directory, and for each BAM file,
    call editing levels at given sites using Query_Editing_Level.pl.

Required arguments:
    --in       Directory containing .bam files (with .bam.bai indexes)
    --out      Output directory where results will be stored
    --ref      Reference genome FASTA file (indexed with samtools faidx)
    --sites    BED file of editing sites (format: chr<TAB>pos<TAB>strand)

Optional:
    --help     Print this help message

Dependencies:
    - samtools (in PATH)
    - Query_Editing_Level.pl (in same directory or PATH)
    - parse_pileup.pl (in same directory as Query_Editing_Level.pl)
EOF
    exit;
}

die "Error: --in, --out, --ref, --sites are required.\n" unless ($in and $out and $ref and $sites);

# Normalize directory paths (add trailing slash if needed)
$in  =~ s/\/$//;
$out =~ s/\/$//;
`mkdir -p $out` unless (-d $out);

# Find all .bam files
opendir(my $dh, $in) or die "Cannot open directory $in: $!";
my @bams = grep { /\.bam$/ } readdir($dh);
closedir($dh);

if (!@bams) {
    die "No .bam files found in $in\n";
}

foreach my $bam (@bams) {
    my $inbam = "$in/$bam";
    my $key = $bam;
    $key =~ s/\.bam$//;
    my $outdir = "$out/$key";
    `mkdir -p $outdir` unless (-d $outdir);
    my $outfile = "$outdir/${key}.all.levels.txt";
    my $cmd = "perl " . dirname(__FILE__) . "/Query_Editing_Level.pl $sites $inbam $outfile";
    print "Running: $cmd\n";
    system($cmd) == 0 or warn "Error running: $cmd\n";
}

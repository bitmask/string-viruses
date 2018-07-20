#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $entities_file = "";
my $pairs_file = "";
my $gzipped = 0;

GetOptions ("entities_file=s" => \$entities_file,
            "pairs_file=s"    => \$pairs_file,
            "gzipped"  => \$gzipped)
           or die("Error in command line arguments\n");

my $root = $ARGV[0];
die "$0 --entities_file=<file> --pairs_file=<file> [--gzipped]" unless $root and $entities_file and $pairs_file;

my %serial_type = ();
my %serial_identifier = ();
open IN, "< $entities_file" or die;
while (<IN>) {
	s/\r?\n//;
	my ($serial, $type, $identifier) = split /\t/;
	$serial_type{$serial} = $type;
	$serial_identifier{$serial} = $identifier;
}
close IN;

if ($gzipped) {
    open IN, "gzip -cd $pairs_file |";
}
else {
    open IN, "cat $pairs_file |";
}

open OUT, "| sort -k1,1n -k3,3n -k5,5nr -T . | gzip -9 > ".$root."_cooccur_pairs.tsv.gz";
while (<IN>) {
	s/\r?\n//;
	my ($serial1, $serial2, $score, undef) = split /\t/;
	next unless exists $serial_type{$serial1} and exists $serial_type{$serial2};
	my $type1;
	my $type2;
	my $identifier1;
	my $identifier2;
	if ($serial_type{$serial1} >= $serial_type{$serial2}) {
		$type1 = $serial_type{$serial1};
		$type2 = $serial_type{$serial2};
		$identifier1 = $serial_identifier{$serial1};
		$identifier2 = $serial_identifier{$serial2};
	}
	else {
		$type1 = $serial_type{$serial2};
		$type2 = $serial_type{$serial1};
		$identifier1 = $serial_identifier{$serial2};
		$identifier2 = $serial_identifier{$serial1};
	}
	print OUT $type1, "\t", $identifier1, "\t", $type2, "\t", $identifier2, "\t", $score, "\n";
}
close IN;
close OUT;

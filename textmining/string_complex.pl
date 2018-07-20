#!/usr/bin/perl -w

use strict;


die unless scalar(@ARGV) == 1;
my $root = $ARGV[0];

my %serial_string = ();
open IN, "< string_entities.tsv";
while (<IN>) {
	s/\r?\n//;
	my ($serial, $type, $identifier) = split /\t/;
	$serial_string{$serial} = $type."\t".$identifier if $type >= -1;
}
close IN;

my $prefix = "PMID"; 
$prefix = "OMIM" if $root eq "omim";
$prefix = "SGD" if $root eq "sgd";

open IN, "cat ".$root.".output |";
open OUT, "| gzip -9 > ".$root."_complex.tsv.gz";
while (<IN>) {
	s/\r?\n//;
	my ($document, $paragraph, $sentence, $start, $end, $name, undef, $serial) = split /\t/;
	next unless exists $serial_string{$serial};
	my $string = $serial_string{$serial};
	print OUT $prefix, ":", $document, "\t", $paragraph, "\t", $sentence, "\t", $start, "\t", $end, "\t", $name, "\t", $string, "\n";
}
close IN;
close OUT;

#!/usr/bin/perl -w

use strict;


die unless scalar(@ARGV) == 1;
my $root = $ARGV[0];

open IN, "gzip -cd ../string_10_5_corpus/all.tsv.gz |";
open OUT, "| gzip -9 > ".$root."_info.tsv.gz";
while (<IN>) {
	s/\r?\n//;
	my ($document, undef, undef, undef, $title, undef) = split /\t/;
	my $url = "";
	$url = "http://omim.org/entry/".$1 if $document =~ /^OMIM:([0-9]+)/;
	$url = "http://www.ncbi.nlm.nih.gov/pubmed/".$1 if $document =~ /^PMID:([0-9]+)/;
	print OUT $document, "\t", $document, "\t", $title, "\t", $url, "\n";
}
close IN;
close OUT;

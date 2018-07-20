#!/usr/bin/perl -w

use strict;
use utf8; 
use Encode; 


die unless scalar(@ARGV) == 1;
my $root = $ARGV[0];

open IN, "gzip -cd ../string_10_5_corpus/all.tsv.gz |";
open OUT, "| gzip -9 > ".$root."_abstract.tsv.gz";
while (<IN>) {
	s/\r?\n//;
	my ($document, undef, undef, undef, $title, $abstract) = split /\t/;
	my $url = "";
	$url = "http://omim.org/entry/".$1 if $document =~ /^OMIM:([0-9]+)/;
	$url = "http://www.ncbi.nlm.nih.gov/pubmed/".$1 if $document =~ /^PMID:([0-9]+)/;
    if (length($abstract) > 2000) {
        $abstract = substr($abstract, 0, 2000);
        $abstract .= "...";
    }
    $abstract =~ s/\\//g;
    $abstract = encode( "UTF-8", $abstract);
	print OUT $document, "\t", $abstract, "\n";
}
close IN;
close OUT;

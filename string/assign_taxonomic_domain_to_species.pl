#! /usr/bin/perl -w

use strict;
use warnings;

my %parent_of_node = ();

open (FH, "$ENV{FREEZE_DIR}/ncbi_taxonomy/nodes.dmp") or die "cannot read file 'nodes.dmp'!\n";

while (<FH>) {
    my ($node, undef, $parent) = split;
    $parent_of_node{$node} = $parent;
}

while (<STDIN>) {

    chomp;
    next if /\A\#/;
    my ($taxon) = split;

    print "$taxon\t";

    my $domain = "UNKNOWN_DOMAIN";

    while (1) {
        last unless exists $parent_of_node{$taxon};
        last if $parent_of_node{$taxon} == $taxon;
        $domain = "bacteria"  if $taxon == 2;
        $domain = "archaea"   if $taxon == 2157;
        $domain = "eukaryota" if $taxon == 2759;
        $domain = "viruses" if $taxon == 10239;
        $taxon = $parent_of_node{$taxon};
    }

    print "$domain\n";
}


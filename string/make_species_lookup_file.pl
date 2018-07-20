#! /usr/bin/perl -w

use strict;
use warnings;

my %valid_species = ();

## first, read in the official species.names.txt file from stdin, to find out which species are
## to be considered.

while (<STDIN>) {
    chomp;
    next if /\A\#/;
    my ($taxon) = split;
    $valid_species{$taxon} = 1;
}

## next, go through all the FASTA files, and connect identifiers to species.

foreach my $species (sort {$a <=> $b} keys %valid_species) {

    my %identifiers_this_species = ();

    my $filename = "$ENV{DERIVED_DIR}/fastafiles/$species.fa";

    open (FH, "< $filename") or die "cannot open '$filename'\n";

    while (<FH>) {

        next unless (/\A\>/);

        if (/\A\>(\d+)\.([A-Za-z0-9\.\-\/_]+)\s/) {
            my $identifier = $2;
            if (exists $identifiers_this_species{$identifier}) {
                die "ERROR: duplicate identifier '$identifier'!\n";
            }
            print "$identifier\t$species\n";
            $identifiers_this_species{$identifier} = 1;

        } else {

            chomp;
            print STDERR "failed to parse '$_'\n";
        }
    }
}


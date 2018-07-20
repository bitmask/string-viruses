#! /usr/bin/perl -w

use strict;
use warnings;

## expect the 'species.lookup' file on stdin.

my %identifiers_per_species = ();

while (<STDIN>) {

    chomp;
    next if /\A\#/;

    my ($identifier, $species) = split;

    if (exists $identifiers_per_species{$species}{$identifier}) {
        warn "WARNING: duplicate protein '$species.$identifier'!\n";
    }

    $identifiers_per_species{$species}{$identifier} = 1;
}

my $current_shorthand = 1;

foreach my $species (sort {$a <=> $b} keys %identifiers_per_species) {

    foreach my $identifier (sort {$a cmp $b} keys %{$identifiers_per_species{$species}}) {

        print "$species.$identifier\t$current_shorthand\n";

        $current_shorthand += 1;
    }
}


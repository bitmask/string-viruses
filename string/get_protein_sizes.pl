#! /usr/bin/perl -w

use strict;
use warnings;

my %identifiers_per_species = ();

open (FH, "species.lookup") or die "cannot read file 'species.lookup'\n";

while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($identifier, $species) = split;
    $identifiers_per_species{$species}{$identifier} = 1;
}

foreach my $species (sort {$a <=> $b} keys %identifiers_per_species) {

    my $fastafile = "$ENV{DERIVED_DIR}/fastafiles/$species.fa";
    open (FH, $fastafile) or die "cannot read file '$fastafile'\n";

    my $current_length = 0;
    my $current_identifier = undef;
    while (<FH>) {
        if (/\A\>/) {
            if (/\A\>(\d+)\.([a-zA-Z0-9_\.\-\/]+)/) {
                print "$species.$current_identifier\t$current_length\n" if defined $current_identifier;
                $current_identifier = $2;
                $current_length = 0;
                die "unknown protein '$species.$current_identifier'\n"
                  unless exists $identifiers_per_species{$species}{$current_identifier};
            } else {
                die "cannot parse '$_'\n";
            }
        } else {
            s/[^A-Za-z]//g;
            $current_length += length;
        }
    }
    print "$species.$current_identifier\t$current_length\n";
}


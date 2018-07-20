#! /usr/bin/perl -w

use strict;
use warnings;

use lib '.';
use lib '../maintenance/';
use Swiss_CRC64;

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

    my %sequences = ();
    my $current_identifier = undef;

    while (<FH>) {
        if (/\A\>/) {
            if (/\A\>([a-zA-Z0-9_\.\-\/]+)/) {
                $current_identifier = $1;
                my ($sp, $id) = split '\.', $current_identifier, 2;
                die "unknown protein '$species.$id'\n" unless exists $identifiers_per_species{$species}{$id};
            } else {
                die "cannot parse '$_'\n";
            }
        } else {
            my $sequence = uc $_;
            $sequence =~ s/[^A-Z]//g;
            $sequences{$current_identifier} .= $sequence;
        }
    }

    foreach my $identifier (sort {$a cmp $b} keys %sequences) {
        my $crc_string = Swiss_CRC64::crc64 ($sequences{$identifier});
        print "$identifier\t$crc_string\n";
    }
}

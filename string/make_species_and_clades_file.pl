#! /usr/bin/perl -w

use strict;
use warnings;

open (FH, "species.names.txt") or die "cannot read file 'species.names.txt'\n";

my %compact_name_of_taxon = ();
my %taxons_per_redundant_name = ();

while (<FH>) {

    chomp;
    next if /\A\#/;

    my ($taxon, $type, $name_official, $name_compact) = split /\t/;

    $compact_name_of_taxon{$taxon} = $name_compact;

    my @words_in_compact_name = split /\s+/, $name_compact;
    my $redundant_name = "$words_in_compact_name[0] $words_in_compact_name[1]";

    ## list additional groupings you may want to have below:

    $redundant_name = "Homo sapiens" if $redundant_name eq "Pan troglodytes";
    $redundant_name = "Mus musculus" if $redundant_name eq "Rattus norvegicus";
    $redundant_name = "Ciona intestinalis" if $redundant_name eq "Ciona savignyi";
    $redundant_name = "Escherichia coli" if $redundant_name =~ /Shigella/;

    $taxons_per_redundant_name{$redundant_name}{$taxon} = 1;
}

my $clade_counter = 0;

open (FH_CLADES, "> clades_and_names.txt") or die "cannot write to file 'clades_and_names.txt'!\n";

foreach my $redundant_name (sort {$a cmp $b} keys %taxons_per_redundant_name) {

    $clade_counter++;

    print "CL$clade_counter ";
    print join " ", sort {$a <=> $b} keys %{$taxons_per_redundant_name{$redundant_name}};
    print "\n";

    print FH_CLADES "CL$clade_counter $redundant_name\n";

    foreach my $taxon (sort {$a <=> $b} keys %{$taxons_per_redundant_name{$redundant_name}}) {
        print "# $taxon $compact_name_of_taxon{$taxon}\n";
    }
}

close FH_CLADES;

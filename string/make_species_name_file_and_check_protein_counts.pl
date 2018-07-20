#! /usr/bin/perl -w

use strict;
use warnings;

my $derived_dir = $ENV{DERIVED_DIR};
my $freeze_dir = $ENV{FREEZE_DIR};

## first, read the names that the import-scripts of Lars have assigned to organisms ...

my ($filename_input_species_names) = $derived_dir."/string_v11.species.tsv";

my %ncbi_name_of_taxon = ();
my %compact_name_of_taxon = ();
my %import_name_of_taxon = ();
my %type_of_taxon = ();

open (FH, $filename_input_species_names) or die "cannot read file '$filename_input_species_names'!\n";
while (<FH>) {

    chomp;
    next if /\A\#/;
    my ($taxon, $ncbi_name, $compact_name, $import_name, $genome_source, $type, $update_taxid) = split /\t/;
    if (exists $ncbi_name_of_taxon{$taxon}) {
        warn "WARNING: double mention of taxon '$taxon' in file '$filename_input_species_names'!\n";
    }
    $ncbi_name_of_taxon{$taxon} = $ncbi_name;
    $compact_name_of_taxon{$taxon} = $compact_name;
    $import_name_of_taxon{$taxon} = $import_name;
    $type_of_taxon{$taxon} = "core species"; 

    $type_of_taxon{$taxon} = "peripheral species" if $type eq 'periphery';
    $type_of_taxon{$taxon} = "adherent species" if $type eq 'adherent';
}

## next, read which FASTA-files we have and how many proteins they contain ...

my $fasta_files = `ls $derived_dir/fastafiles/\*.fa`;
my @fasta_files = split /\s+/, $fasta_files;

my %nr_proteins_per_organism = ();

foreach my $file (@fasta_files) {
    next unless $file =~ /(\d+)\.fa/;
    my $taxon = $1;
    print STDERR "reading fastafile '$taxon.fa'\n";
    open (FH, "$derived_dir/fastafiles/$taxon.fa") or die "cannot read file '$derived_dir/fasta_files/$taxon.fa'!\n";
    my $count = 0;
    my %identifiers_seen = ();
    while (<FH>) {
        next if /\A\#/;
        if (/\A\>/) {
            unless (/\A\>(\d+)\.([\w\.\-]+)\s/) {
                warn "WARNING: cannot parse '_$'\n";
            }
            my $identifier = $2;
            if (exists $identifiers_seen{$identifier}) {
                warn "WARNING: duplicate identifier '$identifier'!\n";
            }
            $identifiers_seen{$identifier} = 1;
            $count++;
        }
    }
    $nr_proteins_per_organism{$taxon} = $count;
}

print "##taxon\ttype\tname_official\tname_compact\tname_NCBI\tname_imported\tnr_of_loci\n";
print "##\n";

foreach my $taxon (sort { $a <=> $b } keys %nr_proteins_per_organism) {
    my $type = $type_of_taxon{$taxon};
    my $ncbi_name = $ncbi_name_of_taxon{$taxon};
    my $compact_name = $compact_name_of_taxon{$taxon};
    if (not defined $ncbi_name) {
        print STDERR $taxon . "\n";
    }
    else {
        print "$taxon\t$type\t$ncbi_name\t$ncbi_name\t$ncbi_name\t$compact_name\t$nr_proteins_per_organism{$taxon}\n";
    }
}

print STDERR "\ndone.\n\nDon't forget to check the result manually, and shorten/simplify the compact names manually if applicable!!\n";


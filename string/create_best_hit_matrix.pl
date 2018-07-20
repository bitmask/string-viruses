#! /usr/bin/perl -w

use strict;
use warnings;
use DBI;

my $param_secondary_hits_cutoff_bitscore = 70;

my $database_name = $ENV{DATABASE_NAME} or die "ERROR: need database name via environment variable !\n";
my $database_host = $ENV{DATABASE_HOST} or die "ERROR: need database host via environment variable !\n";
my $database_port = $ENV{DATABASE_PORT} or die "ERROR: need database port via environment variable !\n";

my %shorthand_of_protein = ();
my %protein_of_shorthand = ();
my %taxon_of_shorthand = ();
my %identifier_of_shorthand = ();

print STDERR "reading file: 'protein.shorthands.txt' ... ";
open (FH, "protein.shorthands.txt") or die "cannot read file 'protein.shorthands.txt'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($protein, $shorthand) = split;
    $shorthand_of_protein{$protein} = $shorthand;
    $protein_of_shorthand{$shorthand} = $protein;
    my ($taxon, $identifier) = $protein =~ /\A(\d+)\.(.*)\z/;
    $taxon_of_shorthand{$shorthand} = $taxon;
    $identifier_of_shorthand{$shorthand} = $identifier;
}
print STDERR "done.\n";

my $count = keys %protein_of_shorthand;
print STDERR "have values in hash $count \n";
my %core_taxa = ();
my %all_taxa = ();

print STDERR "reading file: 'species.names.txt' ... ";
open (FH, "species.names.txt") or die "cannot read file 'species.names.txt'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($taxon, $type) = split /\t/;

    if ($type =~ /core/) {
        $core_taxa{$taxon} = 1;
    }
    $all_taxa{$taxon} = 1;
}
print STDERR "done.\n";

my @sorted_taxa = sort {$a <=> $b} keys %all_taxa;

print "COPY homology.best_hit_per_species FROM stdin;\n";

#my $last_line = <STDIN>;  ## drop the headerline .. 

my $last_line = <STDIN>; chomp $last_line;

foreach my $query_shorthand (sort {$a <=> $b} keys %protein_of_shorthand) {

    print STDERR "." unless $query_shorthand % 100;
    print STDERR " $query_shorthand " unless $query_shorthand % 10000;

    my $protein = $protein_of_shorthand{$query_shorthand};
    my ($home_taxon) = $protein =~ /\A(\d+)\./;
    
    my %best_hits_per_taxon_bitscore = ();
    my %best_hits_per_taxon_shorthand = ();
    my %best_hits_per_taxon_identifier = ();
    my %best_hits_per_taxon_alignment_length = ();
    my %nr_high_scoring_hits_per_taxon = ();

    my $total_best_hit_bitscore = 0;

    while (1) {

        my ($shorthand_a, $shorthand_b, $taxon_b, $bitscore, $start_a, $end_a, $start_b, $end_b, $length_b) = split /\s+/, $last_line;

        last unless $shorthand_a == $query_shorthand;

        $total_best_hit_bitscore = $bitscore if $bitscore > $total_best_hit_bitscore;

        $nr_high_scoring_hits_per_taxon{$taxon_b} += 1 if $bitscore >= $param_secondary_hits_cutoff_bitscore;

        if (exists $best_hits_per_taxon_bitscore{$taxon_b}) {
            unless ($bitscore > $best_hits_per_taxon_bitscore{$taxon_b}) {
                $last_line = <STDIN>; chomp $last_line;
                next;
            }
        }

        my $alignment_length = (($end_a - $start_a) + ($end_b - $start_b)) / 2;
        $alignment_length = int ($alignment_length);
        unless ($alignment_length > 0) {
            print STDERR "WARNING: strange alignment for '$query_shorthand', '$shorthand_b'\n";
            $alignment_length = 0;
        }

        $best_hits_per_taxon_shorthand{$taxon_b} = $shorthand_b;
        $best_hits_per_taxon_bitscore{$taxon_b} = $bitscore;
        $best_hits_per_taxon_identifier{$taxon_b} = $identifier_of_shorthand{$shorthand_b};
        $best_hits_per_taxon_alignment_length{$taxon_b} = $alignment_length;

        $last_line = <STDIN>; 
        chomp $last_line;
    }

    foreach my $taxon (@sorted_taxa) {

        my $valid_comparison = 1;
        $valid_comparison = 0 unless ((exists $core_taxa{$home_taxon}) or (exists $core_taxa{$taxon}));
        $valid_comparison = 1 if $taxon == $home_taxon;

        next unless $valid_comparison;

        my $nr_high_scoring_hits = 0;
        my $best_hit_shorthand = 0;
        my $best_hit_identifier = "-";
        my $best_hit_bitscore = 0;
        my $best_hit_normscore = 0;
        my $best_hit_alignment_length = 0;

        if (exists $best_hits_per_taxon_shorthand{$taxon}) {

            $best_hit_shorthand = $best_hits_per_taxon_shorthand{$taxon};
            $best_hit_identifier = $best_hits_per_taxon_identifier{$taxon};
            $best_hit_bitscore = $best_hits_per_taxon_bitscore{$taxon};
            $best_hit_normscore = $best_hit_bitscore / $total_best_hit_bitscore;
            $best_hit_alignment_length = $best_hits_per_taxon_alignment_length{$taxon};

            ## could be a best hit, but none high-scoring:

            $nr_high_scoring_hits = $nr_high_scoring_hits_per_taxon{$taxon} if exists $nr_high_scoring_hits_per_taxon{$taxon};
        }

        next unless $best_hit_shorthand > 0;    ## optional - use this line to limit the db-size, by not storing non-hits.
                                                ## but, if using this: mind any downstream code that might rely on non-hits
                                                ## being stored as zeros !!

        print "$query_shorthand\t";
        print "$taxon\t";
        print "$nr_high_scoring_hits\t";
        print "$best_hit_shorthand\t";
        print "$best_hit_identifier\t";
        print "$best_hit_bitscore\t";
        print "$best_hit_normscore\t";
        print "$best_hit_alignment_length\n";
    }
}

print "\\.\n";
print STDERR "done.\n";


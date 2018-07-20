#! /usr/bin/perl -w

use strict;
use warnings;
use POSIX;

##############################################################################################
## first, learn about what proteins from which species we are supposed to parse ...
##############################################################################################

my %proteins = ();
my %proteins_per_species = ();
my %species = ();

print STDERR "now processing 'species.lookup' ... ";
open (FH, "< species.lookup") or die "cannot open inputfile 'species.lookup'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($identifier, $species) = split;
    my $protein = "$species.$identifier";
    $proteins{$protein} = 1;
    $species {$species} = 1;
    $proteins_per_species{$species}{$protein} = 1;
}
close FH;
print STDERR "done.\n";

my %shorthand_of_protein = ();

print STDERR "now reading protein-shorthands ...";
open (FH, "protein.shorthands.txt") or die "cannot read file 'protein.shorthands.txt'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($protein, $shorthand) = split;
    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in file 'protein.shorthands.txt'!\n";
    }
    $shorthand_of_protein{$protein} = $shorthand;
}
print STDERR "done.\n";

## now, for each species, parse the FASTA_FILE and get the protein sequence.

print "COPY items.proteins_sequences FROM stdin;\n";

foreach my $species (sort {$a <=> $b} keys %species) {

    my %sequence_of_protein = ();

    my $filename = "$ENV{DERIVED_DIR}/fastafiles/$species.fa";
    open (FH, "< $filename") or die "cannot read file '$filename'\n";

    my $last_protein = undef;

    while (<FH>) {

        chomp;

	if (/\A\>([a-zA-Z0-9_\.\-\/]+)/) {
            my $id = $1;
            $last_protein = "$id";
            unless (exists $proteins{$last_protein}) {
                warn "WARNING: unknown identifier '$id' in species '$species' !\n";
            }
            $sequence_of_protein{$last_protein} = undef;
        } else {
	    if (/\A\>/) {
		warn "WARNING: failed to parse line '$_'!\n";
	    }
            if ($last_protein) {
                my $sequence_bit = uc $_;
                $sequence_bit =~ s/[^A-Z]//g;
                $sequence_of_protein{$last_protein} .= $sequence_bit;
            }
        }

    }

    foreach my $protein (sort { $shorthand_of_protein{$a} <=> $shorthand_of_protein{$b} } keys %sequence_of_protein) {

        my $sequence = $sequence_of_protein{$protein};
        my $shorthand = $shorthand_of_protein{$protein};

        print "$shorthand\t$sequence\n";
    }
}

print "\\.\n";

## that's it.

close STDERR;
close STDOUT;

POSIX::_exit(0);    ## this line: to prevent Perl from going through all its clean-up code.
                    ##            [exits the script the hard way. Destructors etc. will not be called,
                    ##             file handles may not be closed, buffers not flushed. Beware ...]
                    ##
                    ## But: cleaning up does not help anyone, and is darn slow when you have large data structures, so ...


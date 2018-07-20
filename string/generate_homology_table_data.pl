#! /usr/bin/perl -w

use strict;
use warnings;

my %shorthand_of_protein = ();
my %taxon_of_shorthand = ();

my $file_id = shift @ARGV;

my $param_bitscore_cutoff = 50;

my $data_dir = $ENV{DERIVED_DIR};

open (FH, "$data_dir/protein.shorthands.txt") or die "cannot read file 'protein.shorthands.txt'\n";



while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($protein, $shorthand) = split;
    $shorthand_of_protein{$protein} = $shorthand;
    my ($taxon, $identifier) = $protein =~ /\A(\d+)\.(.*)\z/;
    $taxon_of_shorthand{$shorthand} = $taxon;
}

my %core_species = ();
my %all_species = ();

open (FH, "$data_dir/species.names.raw.txt") or die "cannot read file 'species.names.txt'!\n";

while (<FH>) {

    chomp;
    next if /\A\#/;
    my ($species, $type) = split /\t/;

    if ($type =~ /core/) {
        $core_species{$species} = 1;
    }
    $all_species{$species} = 1;
}

my %size_of_protein = ();

open (FH, "$data_dir/protein_sizes.txt") or die "cannot read file 'protein_sizes.txt'\n";

while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($protein, $size) = split;
    if ($size > 32767) {
        warn "WARNING: enforcing a size ceiling: '$_'\n";
        $size = 32767;
    }
    $size_of_protein{$protein} = $size;
}

my %illegal_proteins_skipped = ();
my $already_warned_of_illegal_proteins = 0;

my $file_name = "$ENV{FREEZE_DIR}/SIMAP/$file_id";
my $out_file_name = "$ENV{DERIVED_DIR}/$file_id.parsed_unsorted.tsv.inproc";
my $fin_file_name = "$ENV{DERIVED_DIR}/$file_id.parsed_unsorted.tsv";

unless (-e $file_name) {
    print STDERR "no file with id '$file_id', skipping ...\n";
    exit;
}

open (FH, "gzip -cd $file_name |" ) or die "cannot read file '$file_name'!\n";
open (FH_OUT, "> $out_file_name " ) or die "cannot read file '$out_file_name'!\n";

print STDERR "currently reading file: '$file_name'!\n";

while (<FH>) {
    
    chomp;

    next if /\A\#/;
    
    #my ($protein1, $protein2, $bitscore, $identity, $similarity, $start1, $end1, $start2, $end2, $should_be_empty) = split;
    my ($protein1, $protein2, $unused, $bitscore, $evalue, $identity, $similarity, $start1, $end1, $start2, $end2, $should_be_empty) = split;

    $bitscore =~ s/\,/\./g;

    next unless $bitscore >= $param_bitscore_cutoff;

    $protein1 =~ s/TCOGS2\://g;
    $protein1 =~ s/\:\d+\:pep//g;
    $protein1 =~ s/\:pep//g;
    $protein1 =~ s/\:pseudogenic_transcript//g;
    $protein1 =~ s/PAC\:/PAC_/g;
    $protein1 =~ s/[\(\)]//g;
    $protein1 =~ s/(.)\>/$1/g;

    $protein2 =~ s/TCOGS2\://g;
    $protein2 =~ s/\:\d+\:pep//g;
    $protein2 =~ s/\:pep//g;
    $protein2 =~ s/\:pseudogenic_transcript//g;
    $protein2 =~ s/PAC\:/PAC_/g;
    $protein2 =~ s/[\(\)]//g;
    $protein2 =~ s/(.)\>/$1/g;

    if (($start1 > 32767)
        or ($end1 > 32767)
        or ($start2 > 32767)
        or ($end2 > 32767))
    {
        warn "WARNING: enforcing a ceiling on values beyond 16bits: '$_'\n";
        $start1 = 32767 if $start1 > 32767;
        $end1 = 32767 if $end1 > 32767;
        $start2 = 32767 if $start2 > 32767;
        $end2 = 32767 if $end2 > 32767;
    }
    
    unless ((defined $end2) and (not defined $should_be_empty)) {
        warn "ERROR: failed to parse '$protein1' in file '$file_name'!\n";
	next;
    }
    
    my ($species1) = split (/\./, $protein1);
    my ($species2) = split (/\./, $protein2);
    
    my $line_is_relevant = 0;
    if ($species1 == $species2) {
        $line_is_relevant = 1; 
    } elsif (exists $core_species{$species1}) {
        $line_is_relevant = 1; 
    } elsif (exists $core_species{$species2}) {
        $line_is_relevant = 1;
    }
    
    next unless $line_is_relevant;
    
    my $proteins_are_valid = 1;
    
    my $shorthand1 = $shorthand_of_protein{$protein1};
    my $shorthand2 = $shorthand_of_protein{$protein2};
    
    $proteins_are_valid = 0 unless defined $shorthand1;
    $proteins_are_valid = 0 unless defined $shorthand2;
    
    unless ($proteins_are_valid) {
        unless ($already_warned_of_illegal_proteins) {
            $already_warned_of_illegal_proteins = 1;
            warn "WARNING: no shorthand for protein '$protein1' (further warnings suppressed)\n" unless 
                exists $shorthand_of_protein{$protein1};
            warn "WARNING: no shorthand for protein '$protein2' (further warnings supressed)\n" unless 
                exists $shorthand_of_protein{$protein2};   
        }
        $illegal_proteins_skipped{$protein1} = 1 unless exists $shorthand_of_protein{$protein1};
        $illegal_proteins_skipped{$protein2} = 1 unless exists $shorthand_of_protein{$protein2};
        next;
    }
    
    my $size_protein_1 = $size_of_protein{$protein1};
    my $size_protein_2 = $size_of_protein{$protein2};
   
    my $taxon1 = $taxon_of_shorthand{$shorthand1};
    my $taxon2 = $taxon_of_shorthand{$shorthand2};
 
    ## print both directions here - this will be made nonredundant later !
    
    print FH_OUT "$shorthand1\t$shorthand2\t$taxon2\t$bitscore\t$start1\t$end1\t$start2\t$end2\t$size_protein_2\n";
    print FH_OUT "$shorthand2\t$shorthand1\t$taxon1\t$bitscore\t$start2\t$end2\t$start1\t$end1\t$size_protein_1\n";
}

close FH_OUT;

system("mv $out_file_name $fin_file_name");

my $nr_illegal_proteins_encountered = scalar keys %illegal_proteins_skipped;
print STDERR "$nr_illegal_proteins_encountered illegal/unknown proteins were encountered (and skipped).\n";
if ($nr_illegal_proteins_encountered) {
    print STDERR "These were:\n";
    foreach my $protein (sort {$a cmp $b} keys %illegal_proteins_skipped) {
	print STDERR "$protein\n";
    }
}

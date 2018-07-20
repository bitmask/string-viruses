#! /usr/bin/perl -w

use strict;
use warnings;
use Encode;
use POSIX;
use Sys::Hostname;

## main items;

my %species = ();
my %proteins = ();
my %genes = ();
my %runs = ();
my %orthgroups = ();

## shorthands

my %shorthand_of_protein = ();
my %shorthand_of_orthgroup = ();
my %shorthand_of_gene = ();

## features, connections:

my %species_of_run = ();
my %contig_of_run = ();
my %genes_per_run = ();
my %shorthands_per_run = ();
my %genes_per_protein = ();
my %proteins_per_orthgroup = ();
my %orthgroups_per_protein = ();
my %orthgroups_per_run = ();
my %runs_per_orthgroups = ();
my %preferred_name_of_protein = ();
my %aliases_per_protein = ();
my %proteins_per_alias = ();
my %annotation_of_protein = ();
my %description_lines_per_protein = ();
my %function_line_of_protein = ();
my %annotation_of_orthgroup = ();
my %funccats_per_orthgroup = ();
my %official_name_of_species = ();
my %compact_name_of_species = ();
my %type_of_species = ();
my %kingdom_of_species = ();
my %size_of_protein = ();
my %sequence_of_protein = ();
my %end_position_of_gene = ();
my %orientation_of_gene = ();
my %counts_per_species_per_orthgroup = ();

###############################################################################################
## get the organisms we are to work with
###############################################################################################

print STDERR "learning about the organisms ...";
open (FH, "species.names.txt") or die "cannot read file 'species.names.txt'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($taxon, $type, $name_official, $name_compact) = split /\t/;
    my $type_compact = "periphery";
    $type_compact = "core" if $type =~ /core/;
    $type_compact = "adherent" if $type =~ /adherent/;
    $official_name_of_species{$taxon} = $name_official;
    $compact_name_of_species{$taxon} = $name_compact;
    $type_of_species{$taxon} = $type_compact;
    $species{$taxon} = 1;
}
open (FH, "species.taxonomic_domain.txt") or die "cannot read file 'species.taxonomic_domain.txt'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($taxon, $kingdom) = split;
    unless (exists $species{$taxon}) {
        warn "WARNING: unknown taxon '$taxon' in file 'species.taxonomic_domain.txt'\n";
    }
    $kingdom_of_species{$taxon} = $kingdom;
}

print STDERR "done.\n";

##############################################################################################
## next learn about protein identifiers, and shorthands
##############################################################################################

print STDERR "now processing 'species.lookup' ... ";
open (FH, "species.lookup") or die "cannot open inputfile 'species.lookup'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($identifier, $species) = split;
    my $protein = "$species.$identifier";
    $proteins{$protein} = 1;
    print STDERR " warning: unclear species '$species' in species.lookup!\n" unless exists $species{$species};
}
close FH;
print STDERR "done.\n";

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

###############################################################################################
## now process the orthologous groups.
###############################################################################################

print STDERR "now processing orthologous groups ... ";
open (FH, "final_orthgroups.txt") or die "cannot open inputfile 'final_orthgroups.txt'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($orthgroup, $protein, $start, $end) = split;
    $orthgroups{$orthgroup} = 1;
    if (exists $proteins_per_orthgroup{$orthgroup}{$protein}{$start}) {
        #warn "WARNING: redundant protein/start in 'final_orthgroups.txt': '$_'\n";
        my $old_end = $proteins_per_orthgroup{$orthgroup}{$protein}{$start};
        next unless $end > $old_end;
    }
    $proteins_per_orthgroup{$orthgroup}{$protein}{$start} = $end;
    $orthgroups_per_protein{$protein}{$orthgroup}{$start} = $end;
}
close FH;

print STDERR "done.\n";

###############################################################################################
## next process associated information about orthologous groups.
###############################################################################################

print STDERR "now parsing associated information about orthologous groups ... ";

open (FH, "cat COGs.txt NOGs.txt |") or die "cannot open inputfile 'cat COGs.txt NOGs.txt |'\n";
while (<FH>) {
    chomp;
    next unless (/\A\[([A-Z]+)\]\s+(\w+)\s+(.+)/);
    my $COG = $2;
    my @funccats = split //, $1;
    my $description = $3;
    $description =~ s/[\n\r\']//g;
    if ((length $description) > 990) {
        $description = substr $description, 0, 990;
        $description .= "[...]";
    }
    unless (exists $orthgroups{$COG}) {
        warn "WARNING: unknown orthgroup '$COG' in COGs.txt!\n";
        next;
    }
    $annotation_of_orthgroup{$COG} = $description;
    foreach my $funccat (@funccats) { $funccats_per_orthgroup{$COG}{$funccat} = 1; }
}

print STDERR "done.\n";

#################################################################################################
## checksums ... 
#################################################################################################

my %proteins_per_crc = ();
my %crc_of_protein = ();

print STDERR "now importing protein checksums ...";
open (FH, "proteins.checksums.txt") or die "cannot read file 'proteins.checksums.txt'!\n";
while (<FH>) {
    chomp; next if /\A\#/;
    my ($protein, $checksum) = split;
    $proteins_per_crc{$checksum}{$protein} = 1;
    $crc_of_protein{$protein} = $checksum;
}
print STDERR "done.\n";

#################################################################################################
## now read information about protein sizes.
#################################################################################################

print STDERR "now processing protein sizes ... ";
open (FH, "protein_sizes.txt") or die "cannot open inputfile 'protein_sizes.txt'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($protein, $size) = split;
    $size_of_protein{$protein} = $size;
}
close FH;

print STDERR "done.\n";

#################################################################################################
## now read in information about aliases and preferred names.
#################################################################################################

print STDERR "now processing aliases and preferred names ...";
open (FH, "aliases_and_descriptions/alias_best.tsv") or die "cannot read file 'aliases_and_description/alias_best.tsv'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($taxon, $identifier, $alias, $source) = split /\t/;
    my $protein = "$taxon.$identifier";
    $alias =~ s/\\//g;
    if ($alias =~ /\A\"(.+)\"\z/) { $alias = $1; }

    next if (length($alias) > 48);
    next if (length($alias) < 1);

    if (exists $preferred_name_of_protein{$protein}) {
        #warn "WARNING: redundant protein '$protein' in 'alias_best.tsv'!\n";
    }

    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in 'alias_best.tsv'!\n";
        next;
    }

    $preferred_name_of_protein{$protein} = $alias;
}
open (FH, "aliases_and_descriptions/alias_all.tsv") or die "cannot read file 'aliases_and_description/alias_all.tsv'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($taxon, $identifier, $alias, $sources) = split /\t/;
    my $protein = "$taxon.$identifier";
    $alias =~ s/\\//g;
    if ($alias =~ /\A\"(.+)\"\z/) { $alias = $1; }

    next if (length($alias) > 48);
    next if (length($alias) < 1);


    if (exists $aliases_per_protein{$protein}{$alias}) {
        unless ($alias eq $identifier) {
            #warn "WARNING: redundant alias '$alias - $protein' in 'alias_all.tsv'!\n";
        }
    }
    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in 'alias_all.tsv'!\n";
        next;
    }
    my @sources_split = split /\s+/, $sources;
    foreach my $source (@sources_split) {
        $aliases_per_protein{$protein}{$alias}{$source} = 1;
        $proteins_per_alias{$alias}{$protein}{$source} = 1;
        if ($source =~ /UniProt_DE/) {
            $description_lines_per_protein{$protein}{$alias} = 1;      ## '1' here because low priority: not from the 'text_all' file.
        }
    }
}

## now, for those proteins that did not get a preferred name, and that do come from ensembl, give them the ENSGxxx name as the preferred
## name:

foreach my $protein (keys %shorthand_of_protein) {
    my ($taxon, $identifier) = $protein =~ /\A(\d+)\.(.*)\z/;
    next if exists $preferred_name_of_protein{$protein};
    next unless $identifier =~ /\AENS.*P0\d+\z/;
    foreach my $alias (keys %{$aliases_per_protein{$protein}}) {
	next unless $alias =~ /\AENS.*G0\d+\z/;
	$preferred_name_of_protein{$protein} = $alias;
    }
}

## and, for those proteins that failed to get their own identifier as an alias, add that here:

foreach my $protein (keys %proteins) {

    my ($species, $identifier) = $protein =~ /\A(\d+)\.(.+)\z/;
    next if exists $proteins_per_alias{$identifier}{$protein};
    $aliases_per_protein{$protein}{$identifier}{"direct"} = 1;
    $proteins_per_alias{$identifier}{$protein}{"direct"} = 1;
}

print STDERR "done.\n";

###############################################################################################
## now read in annotations.
###############################################################################################

print STDERR "now importing annotations ...";

open (FH, "aliases_and_descriptions/text_all.tsv") or die "cannot read file 'aliases_and_descriptions/text_all.tsv'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;

    my ($species, $identifier, $annotation, $sources) = split /\t/;
    my $protein = "$species.$identifier";

    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in 'text_all.tsv'!\n";
        next;
    }

    $annotation =~ s/\.\z//;
    my @sources = split /\s/, $sources;
    foreach my $source (@sources) {
        if ($source =~ /UniProt_DE/) {
            $description_lines_per_protein{$protein}{$annotation} = 3;   ## '3' here because this should rank highest ... 
        }
        if ($source eq "Refseq") {
            $description_lines_per_protein{$protein}{$annotation} = 2; 
        }
        if (($source eq "BLAST_UniProt_CC") or ($source eq "Ensembl_UniProt_CC")) {
            my $previous_function_line = "-";
            $previous_function_line = $function_line_of_protein{$protein} if exists $function_line_of_protein{$protein};
            my $length_previous_function_line = length $previous_function_line;
            my $this_length = length $annotation;
            next unless $this_length > $length_previous_function_line;
            $function_line_of_protein{$protein} = $annotation;
        }
    }
}

open (FH, "aliases_and_descriptions/text_best.tsv") or die "cannot read file 'aliases_and_descriptions/text_best.tsv'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    
    my ($species, $identifier, $annotation, $sources) = split /\t/;
    $annotation =~ s/\.\z//;
    
    my $protein = "$species.$identifier";

    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in 'text_best.tsv'!\n";
        next;
    }

    if ($sources =~ /BLAST_UniProt_CC/) {
        if (exists $description_lines_per_protein{$protein}) {
            my @sorted_descriptions = sort {
                my $a_value = $description_lines_per_protein{$protein}{$a};
                my $b_value = $description_lines_per_protein{$protein}{$b};
                return $b_value <=> $a_value unless $b_value == $a_value;
                $a_value = length $a;
                $b_value = length $b;
                return $b_value <=> $a_value; } keys %{$description_lines_per_protein{$protein}};
            my $best_description = shift @sorted_descriptions;
            $annotation = $best_description . "; " . $annotation;
        }
    } else {
        if (exists $function_line_of_protein{$protein}) {
            my $functions = $function_line_of_protein{$protein};
            $annotation = $annotation . "; " . $functions;
        }
    }

    $annotation =~ s/\\//g;
    chop $annotation if $annotation =~ /\.\z/;

    my $encoded_annotation = encode ("utf8", $annotation);

    if (exists $annotation_of_protein{$protein}) {
        warn "WARNING: double mention of protein '$protein' in 'text_best.tsv' ??\n";
    }

    $annotation_of_protein{$protein} = $encoded_annotation;

}
print STDERR "done.\n";

#################################################################################################
## next read in the chromosome positions.
## [note: this is only read to make sure we have all genes. Runs and consolidated chromosome-
##  positions are read further below, from a derived file made by Manuel Stark's scripts.]
#################################################################################################

print STDERR "now reading file gene_position.tsv ...";
open (FH, "gene_position.tsv") or die "cannot read file 'gene_position.tsv'!\n";
while (<FH>) {
    chomp;
    next if /\A\#/;

    my ($species, $identifier, $contig, $orientation, $start, $end) = split /\t/;

    my $protein = "$species.$identifier";

    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in 'gene_position.tsv'!\n";
	next;
    }

    unless (exists $shorthand_of_protein{$protein}) {
        warn "WARNING: unknown protein '$protein' in 'gene_position.tsv'!\n";
	next;
    }

    my $unique_position = $start;
    $unique_position = $end if $end < $start;

    ## 'unique position is need for cases where the same protein-name occurs several times on a contig (identical ORFs, probably .. )

    my $gene = "$protein.$contig.$unique_position";
    if (exists $genes{$gene}) {
        #warn "WARNING: redundant gene '$gene' in 'gene_position.tsv'\n";
    }

    $genes{$gene} = 1;
    $genes_per_protein{$protein}{$gene} = 1;
    $end_position_of_gene{$gene} = $end;
}
print STDERR "done.\n";

##############################################################################################
## now read in the runs
##
##############################################################################################

print STDERR "now reading in the runs ... ";
open (FH, "runs.genes.txt") or die "cannot read file 'runs.genes.txt'\n";
while (<FH>) {
    chomp;
    next if /\A\#/;
    my ($run, $gene_in_manuel_notation, $orientation) = split;
    my ($species, $identifier, $chromosome, $start_position, $end_position) = split /\#/, $gene_in_manuel_notation;
    unless ((defined $start_position) and (defined $orientation) and (defined $end_position)) {
        warn "WARNING: unable to parse '$_' in file 'runs.genes.txt'\n";
	next;
    }
    my $protein = "$species.$identifier";
    unless (exists $proteins{$protein}) {
        warn "WARNING: unknown protein '$protein' in file 'runs.genes.txt'\n";
	next;
    }
    unless (exists $shorthand_of_protein{$protein}) {
        warn "ERROR: no shorthand for protein '$protein' in file 'runs.genes.txt'!\n";
	next;
    }
    my $gene = "$species.$identifier.$chromosome.$start_position";

    unless (exists $genes{$gene}) {
	warn "WARNING: unknown gene '$gene' in file 'runs.genes.txt'!\n";
	next;
    }

    $genes{$gene} = 1;
    $genes_per_protein{$protein}{$gene} = 1;
    $end_position_of_gene{$gene} = $end_position;      ## we just overwrite here, in case something was stored already. 
                                                       ## this is done to make sure that positions and runs are in synch...
    $orientation_of_gene{$gene} = $orientation;

    my $shorthand = $shorthand_of_protein{$protein};
    $species_of_run{$run} = $species;
    $contig_of_run{$run} = $chromosome;
    $genes_per_run{$run}{$gene} = 1;
    $shorthands_per_run{$run}{$shorthand} = 1;
    if (exists $orthgroups_per_protein{$protein}) {
        foreach my $orthgroup (keys %{ $orthgroups_per_protein{$protein}}) {
            $orthgroups_per_run{$run}{$orthgroup} = 1;
        }
    }
    $runs{$run} = 1;
}
print STDERR "done.\n";

## now go through the complete set of proteins and check whether there are any
## that have not yet been assigned to a chromosomal position (These should be largely
## higher eukaryotes). Create genes identifiers for them, and sign them up to
## the list of genes.

foreach my $protein (keys %proteins) {
    next if (exists $genes_per_protein{$protein});
    my $gene = "$protein.unknown.1";         ## we create a dummy gene entry here for cases
    $genes{$gene} = 1;                       ## where there is no chromosome information
    $genes_per_protein{$protein}{$gene} = 1;
}

## connect orthgroups and species ...

foreach my $orthgroup (keys %proteins_per_orthgroup) {
    foreach my $protein (keys %{$proteins_per_orthgroup{$orthgroup}}) {
	my ($species) = $protein =~ /\A(\d+)\./;
	$counts_per_species_per_orthgroup{$orthgroup}{$species} += 1;
    }
}

## now assign internal shorthand-identifiers for orthgroups and genes:

my $gene_shorthand = 1;

foreach my $gene (sort {$a cmp $b} keys %genes) {
    $shorthand_of_gene{$gene} = $gene_shorthand;
    $gene_shorthand += 1;
}

## WATCH !
## the offset below needs to be larger than the number of proteins - this should not be hardcoded here in the long run !

my $orthgroup_shorthand = 100000000;

foreach my $orthgroup (sort {$a cmp $b} keys %orthgroups) {
    $shorthand_of_orthgroup{$orthgroup} = $orthgroup_shorthand;
    $orthgroup_shorthand += 1;
}


print STDERR "now printing to stdout.\n";

## now simply output the stuff by writing out the SQL commands:

## first the main tables ...

print "COPY items.species FROM stdin;\n";
foreach my $species (sort {$a <=> $b} keys %species) {
    my $official_name = $species;
    $official_name = $official_name_of_species{$species} if exists $official_name_of_species{$species};
    my $compact_name = $species;
    $compact_name = $compact_name_of_species{$species} if exists $compact_name_of_species{$species};
    my $type = "unknown";
    $type = $type_of_species{$species} if exists $type_of_species{$species};
    my $kingdom = "unknown";
    $kingdom = $kingdom_of_species{$species} if exists $kingdom_of_species{$species};
    print "$species\t" . "$official_name\t" . "$compact_name\t" . "$kingdom\t" . "$type\n";
}
print "\\.\n";

print "COPY items.genes FROM stdin;\n";
foreach my $gene (sort {$a cmp $b} keys %genes) {
    my ($taxon, $identifier, $contig, $start_position) = $gene =~ /\A(\d+)\.(.+)\.([\w\-]+)\.(\d+)/;
    unless (defined $start_position) {
        warn "WARNING: failed to parse gene '$gene' in \#230958\n";
        next;
    }
    my $protein = "$taxon.$identifier";
    my $size = 10;
    unless (exists $size_of_protein{$protein}) {
        warn "WARNING: no size for protein '$protein'!\n";
    }
    $size = $size_of_protein{$protein} if exists $size_of_protein{$protein};
    my $end_position = 2;
    $end_position = $end_position_of_gene{$gene} if exists $end_position_of_gene{$gene};

    print "$shorthand_of_gene{$gene}\t" . "$gene\t" . "$start_position\t" . "$end_position\t" . "$size\n";

    ## note: the above line does not contain direct information about gene orientation. 
    ##       instead, this is encoded in the positions: the orientation is 'reverse', if the 
    ##       end position is larger than the start position. 
}
print "\\.\n";

print "COPY items.proteins FROM stdin;\n";
my $warn_counter_proteins1 = 0;
my $warn_counter_proteins2 = 0;
foreach my $protein (sort {$shorthand_of_protein{$a} <=> $shorthand_of_protein{$b}} keys %proteins) {
    my ($taxon, $identifier) = $protein =~ /\A(\d+)\.(.+)\z/;
    my $size = 10;
    unless (exists $size_of_protein{$protein}) {
        warn "WARNING: no size for protein '$protein'!\n";
    }
    $size = $size_of_protein{$protein} if exists $size_of_protein{$protein};
    unless (exists $shorthand_of_protein{$protein}) {
        warn "WARNING: no shorthand for protein '$protein'!\n";
    }
    my $shorthand = 0;
    $shorthand = $shorthand_of_protein{$protein} if exists $shorthand_of_protein{$protein};
    unless (exists $crc_of_protein{$protein}) {
        $warn_counter_proteins1 += 1;
        warn "WARNING: no CRC checksum for protein '$protein'!\n" unless $warn_counter_proteins1 > 100;
        warn "WARNING: stopped reporting - too many warnings.\n" if $warn_counter_proteins1 == 101;
    }
    my $crc_checksum = 0;
    $crc_checksum = $crc_of_protein{$protein} if exists $crc_of_protein{$protein};
    my $annotation = "annotation not available";
    $annotation = $annotation_of_protein{$protein} if exists $annotation_of_protein{$protein};
    if (length($annotation) > 598) {
        $annotation = substr $annotation, 0, 592;
        $annotation .= " [...] ";
    }
    unless (exists $preferred_name_of_protein{$protein}) {
        $warn_counter_proteins2 += 1;
        warn "WARNING: no preferred name for protein '$protein'!\n" unless $warn_counter_proteins2 > 100;
        warn "WARNING: stopped reporting - too many warnings.\n" if $warn_counter_proteins2 == 101;
    }
    my $preferred_name = $identifier;
    $preferred_name = $preferred_name_of_protein{$protein} if exists $preferred_name_of_protein{$protein};
    print "$shorthand\t" . "$protein\t" . "$taxon\t" . "$crc_checksum\t" . "$size\t" . "$annotation\t" . "$preferred_name\n";
}
print "\\.\n";

print "COPY items.runs FROM stdin;\n";
foreach my $run (sort {$a <=> $b} keys %runs) {
    die "ERROR: no species for run '$run'!\n" unless exists $species_of_run{$run};
    die "ERROR: no contig for run '$run'!\n"  unless exists $contig_of_run{$run};
    my $species = $species_of_run{$run};
    my $contig  = $contig_of_run{$run};
    print "$run\t" .
	"$species\t" .
	"$contig\n";
}
print "\\.\n";

print "COPY items.orthgroups FROM stdin;\n";
foreach my $orthgroup (sort {$a cmp $b} keys %orthgroups) {
    my $annotation = "non supervised orthologous group";
    $annotation = $annotation_of_orthgroup{$orthgroup} if (exists $annotation_of_orthgroup{$orthgroup});
    my $protein_count = scalar keys %{$proteins_per_orthgroup{$orthgroup}};
    my $species_count = scalar keys %{$counts_per_species_per_orthgroup{$orthgroup}};
    print "$shorthand_of_orthgroup{$orthgroup}\t" . 
	"$orthgroup\t" . 
	"$annotation\t" .
	"$protein_count\t" .
	"$species_count\n";
}
print "\\.\n";

## ... then the connecting tables.

print "COPY items.orthgroups_funccats FROM stdin;\n";
foreach my $orthgroup (sort {$a cmp $b} keys %funccats_per_orthgroup) {
    foreach my $funccat (sort {$a cmp $b} keys %{$funccats_per_orthgroup{$orthgroup}}) {
        print "$shorthand_of_orthgroup{$orthgroup}\t" . 
	    "$funccat\n";
    }
}
print "\\.\n";

print "COPY items.proteins_names FROM stdin;\n";
foreach my $alias (sort {$a cmp $b} keys %proteins_per_alias) {
    foreach my $protein (sort {$a cmp $b} keys %{$proteins_per_alias{$alias}}) { 
        my ($species) = $protein =~ /\A(\d+)\./;
        foreach my $source (sort {$a cmp $b} keys %{$proteins_per_alias{$alias}{$protein}}) {
            my $preferred_name_flag = 'false';
            if (exists $preferred_name_of_protein{$protein}) {
                $preferred_name_flag = 'true' if $preferred_name_of_protein{$protein} eq $alias;
            }
            print "$alias\t" .
                "$shorthand_of_protein{$protein}\t" .
                "$species\t" .
                "$source\t" .
                "$preferred_name_flag\n";
        }
    }
}
print "\\.\n";

print "COPY items.proteins_orthgroups FROM stdin;\n";
foreach my $orthgroup (sort {$a cmp $b} keys %proteins_per_orthgroup) {
    foreach my $protein (sort {$a cmp $b} keys %{$proteins_per_orthgroup{$orthgroup}}) {
	my ($species, $identifier) = $protein =~ /\A(\d+)\.(.+)\z/;
	my $preferred_name = $identifier;
	$preferred_name = $preferred_name_of_protein{$protein} if exists $preferred_name_of_protein{$protein};
        my $annotation = "annotation not available";
        $annotation = $annotation_of_protein{$protein} if exists $annotation_of_protein{$protein};
        if (length($annotation) > 98) {
            $annotation = substr $annotation, 0, 92;
            $annotation .= " [...] ";
        }
        foreach my $startpos (keys %{$proteins_per_orthgroup{$orthgroup}{$protein}}) {
            my $endpos = $proteins_per_orthgroup{$orthgroup}{$protein}{$startpos};
            print "$shorthand_of_orthgroup{$orthgroup}\t" . 
		"$shorthand_of_protein{$protein}\t" . 
                "$protein\t" .
		"$species\t" .
		"$startpos\t" . 
		"$endpos\t" . 
		"$preferred_name\t" .
		"$annotation\t" .
		"\\N\n";
        }
    }
}
print "\\.\n";

print "COPY items.genes_proteins FROM stdin;\n";
foreach my $protein (sort {$shorthand_of_protein{$a} <=> $shorthand_of_protein{$b}} keys %genes_per_protein) {
    foreach my $gene (sort {$shorthand_of_gene{$a} <=> $shorthand_of_gene{$b}} keys %{$genes_per_protein{$protein}}) {
        print "$shorthand_of_protein{$protein}\t";
        print "$shorthand_of_gene{$gene}\n";
    }
}
print "\\.\n";

print "COPY items.runs_genes_proteins FROM stdin;\n";
foreach my $run (sort {$a <=> $b} keys %runs) {
    unless (exists $genes_per_run{$run}) {
        die "ERROR: no gene in run '$run'!\n";
    }
    foreach my $gene (sort {$shorthand_of_gene{$a} <=> $shorthand_of_gene{$b}} keys %{$genes_per_run{$run}}) {

        my ($taxon, $identifier, $contig, $start_position) = $gene =~ /\A(\d+)\.(.+)\.([\w\-]+)\.(\d+)/;
        unless (defined $start_position) {
            warn "WARNING: failed to parse gene '$gene' in \#230959\n";
            next;
        }

        my $end_position = 2;
        $end_position = $end_position_of_gene{$gene} if exists $end_position_of_gene{$gene};
        my $protein = "$taxon.$identifier";

	my $preferred_name = $identifier;
	$preferred_name = $preferred_name_of_protein{$protein} if exists $preferred_name_of_protein{$protein};

        my $annotation = "annotation not available";
        $annotation = $annotation_of_protein{$protein} if exists $annotation_of_protein{$protein};

        if (length($annotation) > 98) {
            $annotation = substr $annotation, 0, 92;
            $annotation .= " [...] ";
        }
        ## by convention: turn positions around if encoded on antisense strand ...
        if (exists $orientation_of_gene{$gene}) {
            if ($orientation_of_gene{$gene} eq "-") { 
                my $dummy = $start_position; $start_position = $end_position; $end_position = $dummy; 
            }
        }
        print "$run\t" .
            "$shorthand_of_gene{$gene}\t" .
            "$shorthand_of_protein{$protein}\t" .
            "$start_position\t" .
            "$end_position\t" .
            "$preferred_name\t" .
	    "$annotation\n";
    }
}
print "\\.\n";

print "COPY items.runs_orthgroups FROM stdin;\n";
foreach my $run (sort {$a <=> $b} keys %orthgroups_per_run) {
    foreach my $orthgroup (sort {$shorthand_of_orthgroup{$a} <=> $shorthand_of_orthgroup{$b}} keys %{$orthgroups_per_run{$run}}) {
        print "$run\t" . "$shorthand_of_orthgroup{$orthgroup}\n";
    }
}
print "\\.\n";

print "COPY items.orthgroups_species FROM stdin;\n";
foreach my $orthgroup (sort {$shorthand_of_orthgroup{$a} <=> $shorthand_of_orthgroup{$b}} keys %counts_per_species_per_orthgroup) {
    foreach my $species (sort {$a <=> $b} keys %species) {
	my $count = 0;
	$count = $counts_per_species_per_orthgroup{$orthgroup}{$species} if 
	    exists $counts_per_species_per_orthgroup{$orthgroup}{$species};
	next if $count < 1;    ## needed for space efficiency at the database
	print "$shorthand_of_orthgroup{$orthgroup}\t" .
	    "$species\t" .
	    "$count\n";
    }
}
print "\\.\n";

print STDERR "\nall done.\n";

## that's it.

close STDERR;
close STDOUT;

POSIX::_exit(0);    ## this line: to prevent Perl from going through all its clean-up code.
                    ##            [exits the script the hard way. Destructors etc. will not be called,
                    ##             file handles may not be closed, buffers not flushed. Beware ...]
                    ##
                    ## But: cleaning up does not help anyone, and is darn slow when you have large data structures, so ...


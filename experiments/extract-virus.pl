
use warnings;
use strict;
use Getopt::Long;

# get virus and cross species interactions, sort them and write to output, zipped
# used to write to individual species files -- but is this necessary?

my $file = "";
my $viruses = "";
my $output_prefix = "";
my $mode = "";
my $nohost = 0;
my $transferformat = 0;
my $ps_file = "";

GetOptions("scores=s" => \$file,
           "viruses=s" => \$viruses,
           "output=s" => \$output_prefix,
           "mode=s" => \$mode,
           "nohost" => \$nohost, # print virus-virus interactions only
           "transferformat" => \$transferformat, # sort viruses and hosts into their own columns, requires proteinshorthands
           "proteinshorthands=s" => \$ps_file, # stores id for each protein
       ) or die "error in command line arguments";

die "$0 --scores=<file> --viruses=<file> --output=<file> --mode=<type>" unless $file and $viruses and $output_prefix and $mode;

my %virus;
open VIR, "$viruses" or die "can't open $viruses for reading $!";
while (<VIR>) {
    my ($taxid, undef) = split /\t/;
    $virus{$taxid} = 1;
}
close VIR;

my %shorthands;
if ($transferformat) {
    print STDERR "opening shorthands...";
    open PS, "$ps_file" or die "can't open $ps_file for reading $!";
    while (<PS>) {
        chomp;
        my ($string, $short) = split /\t/;
        $shorthands{$string} = $short;
    }
    close PS;
    print STDERR "done\n";
}

my $number = 0;
$number = 8 if $mode eq 'experiments';
$number = 12 if $mode eq 'textmining';
$number = 12 if $mode eq 'actions';

my %output;
open FH, "zcat $file |" or die "can't open $file for reading";
while (<FH>) {
    chomp;
    my $line = $_;
    if ($mode eq 'experiments' or $mode eq 'textmining') {
        my ($label, $taxid, $prot1, $prot2, $prob, $pmid, $evidence) = split /\t/;
        if (not defined $virus{$taxid}) {
            if ($taxid =~ /-/) { # then it is cross species, and we want to consider x-y and y-x to be the same id
                my ($taxid1, $taxid2) = split /-/, $taxid;
                if ($transferformat) {
                    my $short1 = $shorthands{$taxid1 . "." . $prot1};
                    warn "no shorthand for $taxid1 $prot1" if not defined $short1;
                    my $short2 = $shorthands{$taxid2 . "." . $prot2};
                    warn "no shorthand for $taxid2 $prot2" if not defined $short2;
                    if (defined $virus{$taxid1}) {
                        push @{$output{$taxid}}, "$taxid1-$taxid2\t$short1\t$short2\t$number\t$number\t$prob\n" unless $nohost;
                    }
                    elsif (defined $virus{$taxid2}) {
                        push @{$output{$taxid}}, "$taxid2-$taxid1\t$short2\t$short1\t$number\t$number\t$prob\n" unless $nohost;
                    }
                    else {
                        die "this should not happen"
                    }
                }
                else {
                    push @{$output{$taxid}}, "$taxid1.$prot1\t$taxid2.$prot2\t$prob\n" unless $nohost;
                }
            }
            # ignore anything that is host-host
        }
        else {
            if ($transferformat) {
                my $continue = 1;
                if (not defined $shorthands{$taxid.".".$prot1}) {
                    print STDERR "$taxid.$prot1 has no shorthand\n";
                    $continue = 0;
                }
                if (not defined $shorthands{$taxid.".".$prot2}) {
                    print STDERR "$taxid.$prot2 has no shorthand\n";
                    $continue = 0;
                }
                if ($continue) {
                    my $short1 = $shorthands{$taxid.".".$prot1};
                    my $short2 = $shorthands{$taxid.".".$prot2};
                    push @{$output{$taxid}}, "$taxid\t$short1\t$short2\t$number\t$number\t$prob\n";
                }
            }
            else {
                push @{$output{$taxid}}, "$taxid.$prot1\t$taxid.$prot2\t$prob\n";
            }
        }
    }
    elsif ($mode eq 'actions') {
        # Replace this with translation script in maintenance????
        my ($label, $taxid, $prot1, $prot2, $score, $pmid, $name1, $name2, $text, $start, $end) = split /\t/;
        if (not defined $virus{$taxid}) {
            if ($taxid =~ /(.*)-(.*)/) {
                # cross species
                my ($taxid1, $taxid2) = ($1, $2);
                    if ($transferformat) {
                        my $short1 = $shorthands{$taxid1.".".$prot1};
                        my $short2 = $shorthands{$taxid2.".".$prot2};
                        # skip NOG and KOG entries, which aren't in shorthands
                        push @{$output{$taxid}}, "$taxid1\t$short1\t$short2\t$label\t$number\t$score\n" if defined $short2;
                    }
            }
        }
        else {
            # virus-virus
            if ($transferformat) {
                my $short1 = $shorthands{$taxid.".".$prot1};
                my $short2 = $shorthands{$taxid.".".$prot2};
                push @{$output{$taxid}}, "$taxid\t$short1\t$short2\t$label\t$number\t$score\n";
            }
        }
    }
}

my $outfile = $output_prefix; # . $taxid;
open OUT, ">$outfile";
foreach my $taxid (keys %output) {
    foreach my $line (@{$output{$taxid}}) {
        print OUT $line;
    }
}

#!/usr/bin/perl -w

use strict;

use Getopt::Long;

my $textmining = "";
my $nlp = "";
my $experiments = "";
my $transfer = "";
my $output = "";  # in the format lars wants for the cytoscape app, or in the format needed to import to the database
my $shorthands = "";
my $prior = 0.041;

GetOptions ("textmining=s" => \$textmining,
            "nlp=s" => \$nlp,
            "experiments=s"    => \$experiments,
            "transfer=s" => \$transfer, # all transfer evidence concatenated together
            "output=s" => \$output,
            "shorthands=s" => \$shorthands,
           ) or die("Error in command line arguments\n");

die "$0 --textmining=<file> --experiments=<file> --nlp=<file> --transfer=<file> --output=<type>" unless $textmining and $experiments and $nlp;

open (my $tm, "<", $textmining) or die "can't open $textmining $!";
open (my $nl, "<", $nlp) or die "can't open $nlp $!";
open (my $ex, "<", $experiments) or die "can't open $experiments $!";
open (my $tr, "<", $transfer) or die "can't open $transfer $!";

my %protein_to_shorthand;
my %shorthand_to_protein; # transfer output is on shorthand
open(my $fh, "<", $shorthands) or die "Can't open shorthands for reading: $!";
while (my $line = <$fh>) {
    chomp $line;
    my ($shorthand, $string) = split "\t", $line;
    $protein_to_shorthand{$string} = $shorthand;
    $shorthand_to_protein{$shorthand} = $string;
}

if ($output eq "db") {
    print "COPY network.node_node_links FROM stdin;\n";
}

my $scores;
$scores = read_scores($tm, "tm", $scores);
$scores = read_nlp_scores($nl, "nl", $scores);
$scores = read_scores($ex, "ex", $scores);
$scores = read_transfer_scores($tr, "tr", $scores);
close $tm;
close $ex;


use Data::Dumper;

sub subtract_prior {
    my ($score) = @_;
    $score = $score - $prior;
    $score = 0 if $score < 0;
    $score = $score / (1-$prior);
    $score = 1 if $score > 1;
    return $score;
}

sub combine_scores {
    my ($scores) = @_;
    #subtract priors
    my @without_priors = ();
    foreach my $s (@$scores) {
        push @without_priors, subtract_prior($s);
    }

    my $mult = 1;
    foreach my $s (@without_priors) {
        $mult *= (1-$s);
    }
    my $final_score = (1 - $mult) * (1-$prior) + $prior;
    $final_score = 1 if $final_score > 1;
    $final_score = 0 if $final_score < 0;
    return $final_score;
}

foreach my $string1 (keys %$scores) {
    foreach my $string2 (keys %{$scores->{$string1}}) {
        if (defined $scores->{$string1}->{$string2}) {
            next if $string1 eq $string2;
            my $ex = 0;
            my $tm = 0;
            my $nl = 0;
            my $te = 0;
            my $tt = 0;
            my $tn = 0;

            my @scores = ();
            $ex = $scores->{$string1}->{$string2}->{"ex"} if defined $scores->{$string1}->{$string2}->{"ex"};
            $tm = $scores->{$string1}->{$string2}->{"tm"} if defined $scores->{$string1}->{$string2}->{"tm"};
            $nl = $scores->{$string1}->{$string2}->{"nl"} if defined $scores->{$string1}->{$string2}->{"nl"};
            #$tm = $nl if $nl > $tm; # tm score is greater of nlp and cooccurrence --- don't do this
            $tm = combine_scores([$tm, $nl]) if $tm > 0 and $nl > 0;
            $tm = $nl if $tm == 0;

            #transfer
            $te = $scores->{$string1}->{$string2}->{"tr8 8"} if defined $scores->{$string1}->{$string2}->{"tr8 8"};
            $tt = $scores->{$string1}->{$string2}->{"tr12 12"} if defined $scores->{$string1}->{$string2}->{"tr12 12"};

            # take largest nlp transfer score
            my $largest = 0;
            foreach my $s (keys(%{$scores->{$string1}->{$string2}})) {
                my $t = $scores->{$string1}->{$string2}->{$s} if $s =~ /tr[A-z]/; 
                $largest = $t if defined $t and $t > $largest;
            }
            $tn = $largest;
            #$tt = $tn if $tn > $tt;
            $tt = combine_scores([$tt, $tn]) if $tt > 0 and $tn > 0;
            $tt = $tn if $tt == 0;

            push @scores, $ex if $ex > 0;
            push @scores, $tm if $tm > 0;
            push @scores, $te if $ex > 0;
            push @scores, $tt if $tt > 0;
            
            my $combined_score = combine_scores(\@scores);
            my $comb = sprintf("%.3f", $combined_score);
            next if $combined_score*1000 < 150;
            if ($output eq "lars") {
                next if $string1 =~ /NOG|COG|KOG/ or $string2 =~ /NOG|COG|KOG/;
                # combine the scores with their transfer scores for Lars
                $ex = combine_scores([$ex, $te]) if $ex > 0 and $te > 0;
                $tm = combine_scores([$tm, $tt]) if $tm > 0 and $tt > 0;
                $ex = $te if $ex == 0;
                $tm = $tt if $tm == 0;
                $ex = sprintf("%.3f", $ex);
                $tm = sprintf("%.3f", $tm);
                print $string1 . "\t" . $string2 . "\t" . 0 . "\t" . 0 . "\t" . 0 . "\t" . 0 . "\t" . 0 . "\t" . $ex*1000 . "\t" . 0 . "\t" . $tm*1000 . "\t" . $comb*1000 . "\n";
                # protein1 protein2 sscore neighborhood fusion cooccurence coexpression experimental database textmining combined_score
            }
            else {
                # keep these scores separate from the transfered scores for the database
                my ($taxid, $protein_id) = split '\.', $string2, 2;
                my $sh1 = "";
                $sh1 = $protein_to_shorthand{$string1} if defined $protein_to_shorthand{$string1};
                my $sh2 = "";
                $sh2 = $protein_to_shorthand{$string2} if defined $protein_to_shorthand{$string2};
                my @ar;
                $ex = sprintf("%.3f", $ex);
                $tm = sprintf("%.3f", $tm);
                $te = sprintf("%.3f", $te);
                $tt = sprintf("%.3f", $tt);
                push @ar, "{8, " . $ex*1000 . "}" if $ex > 0;
                push @ar, "{9, " . $te*1000 . "}" if $te > 0;
                push @ar, "{12, " . $tm*1000 . "}" if $tm > 0;
                push @ar, "{13, " . $tt*1000 . "}" if $tt > 0;
                my $array = "{";
                $array .= join ",", @ar;
                $array .= "}";
                if ($output eq "r") {
                    my ($taxid1, $protein_id1) = split '\.', $string1, 2;
                    print $taxid1 . "\t" . $taxid . "\tex\t" . $ex . "\n" if $ex;
                    print $taxid1 . "\t" . $taxid . "\ttm\t" . $tm . "\n" if $tm;
                    print $taxid1 . "\t" . $taxid . "\tte\t" . $te . "\n" if $te;
                    print $taxid1 . "\t" . $taxid . "\ttt\t" . $tt . "\n" if $tt;
                }
                else {
                    print $sh1 . "\t" . $taxid . "\t" . $sh2 . "\t".$comb * 1000 . "\t" . $array . "\n" if $sh1 and $sh2;
                }
            }
        }
    }
}

sub read_scores {
    my ($fh, $label, $scores) = @_;
    while (<$fh>) {
        chomp;
        my ($string1, $string2, $score) = split "\t", $_;
        $scores->{$string1}->{$string2}->{$label} = $score;

    }
    return $scores;
}

sub read_nlp_scores {
    my ($fh, $label, $scores) = @_;
    while (<$fh>) {
        chomp;
        my ($source, $taxids, $protein1, $protein2, $score) = split "\t", $_;
        
        my ($taxid1, $taxid2) = split "-", $taxids;
        $taxid2 = $taxid1 unless $taxid2;
        my $string1 = "$taxid1.$protein1";
        my $string2 = "$taxid2.$protein2";
        $scores->{$string1}->{$string2}->{$label} = $score;
        $scores->{$string2}->{$string1}->{$label} = $score;  # nlp file is not already symmetric
    }
    return $scores;
}

sub read_transfer_scores {
    my ($fh, $tr, $scores) = @_;
    while (<$fh>) {
        chomp;
        my ($htaxid, $sh1, $sh2, $label, $score, $vtaxid) = split "\t", $_;
        my $string1 = $shorthand_to_protein{$sh1};
        my $string2 = $shorthand_to_protein{$sh2};
        $scores->{$string1}->{$string2}->{$tr.$label} = $score;
        $scores->{$string2}->{$string1}->{$tr.$label} = $score;  # transfer is symmetric?????
    }
    return $scores;
}

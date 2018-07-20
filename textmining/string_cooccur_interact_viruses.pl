#!/usr/bin/perl -w

use strict;


die unless scalar(@ARGV) == 1;
my $root = $ARGV[0];

my %type_identifier_kegg = ();
open IN, "< benchmark_kegg.tsv" or die;
while (<IN>) {
	s/\r?\n//;
	my ($type, $kegg, undef, @identifiers) = split /[\t ]/;
	foreach my $identifier (@identifiers) {
		if ($identifier =~ /^CID[ms][0-9]+$/) {
			$type_identifier_kegg{"-1"}{$identifier}{$kegg} = 1;
		}
		else {
			$type_identifier_kegg{$type}{$identifier}{$kegg} = 1;
		}
	}
}
close IN;

my %type_type_score_prob = ();
my @bincount = ();
my @binscoresum = ();
my @binstatesum = ();
my $type1 = "";
my $type2 = "";

sub monotonize
{
	my $continue = 1;
	while ($continue) {
		$continue = 0;
		for (my $i = 0; $i < $#bincount; $i++) {
			next unless exists $bincount[$i];
			my $j = $i+1;
			$j++ until exists $bincount[$j] or $j > $#bincount;
			last unless exists $bincount[$j];
			if ($binstatesum[$i]/$bincount[$i] <= $binstatesum[$j]/$bincount[$j]) {
				$bincount[$i] += $bincount[$j];
				$binscoresum[$i] += $binscoresum[$j];
				$binstatesum[$i] += $binstatesum[$j];
				delete $bincount[$j];
				delete $binscoresum[$j];
				delete $binstatesum[$j];
				$continue = 1;
			}
		}
	}
	for (my $i = 0; $i <= $#bincount; $i++) {
		$type_type_score_prob{$type1}{$type2}{$binscoresum[$i]/$bincount[$i]} = $binstatesum[$i]/$bincount[$i] if exists $bincount[$i];
	}
}

#open IN, "gzip -cd ".$root."_cooccur_pairs.tsv.gz |";
open IN, "human"; # read the human entries to build the benchmarking curve
my $count = ();
my $scoresum = 0;
my $statesum = 0;
my @scores = ();
my @states = ();
while (<IN>) {
	s/\r?\n//;
	my ($newtype1, $identifier1, $newtype2, $identifier2, $score) = split /\t/;
	if ($newtype1 ne $type1 or $newtype2 ne $type2) {
		monotonize() if $type1 ne "" and $type2 ne "";
		$count = 0;
		$scoresum = 0;
		$statesum = 0;
		@scores = ();
		@states = ();
		@bincount = ();
		@binscoresum = ();
		@binstatesum = ();
		$type1 = $newtype1;
		$type2 = $newtype2;
	}
	
	next unless exists $type_identifier_kegg{$type1} and exists $type_identifier_kegg{$type1}{$identifier1};
	next unless exists $type_identifier_kegg{$type2} and exists $type_identifier_kegg{$type2}{$identifier2};
	my $state = 0;
	foreach my $kegg (keys %{$type_identifier_kegg{$type1}{$identifier1}}) {
		if (exists $type_identifier_kegg{$type2}{$identifier2}{$kegg}) {
			$state = 1;
			last;
		}
	}
	$count++;
	$scoresum += $score;
	$statesum += $state;
	push @states, $state;
	push @scores, $score;
	if ($count >= 100) {
		push @bincount, $count;
		push @binscoresum, $scoresum;
		push @binstatesum, $statesum;
		$count--;
		$scoresum -= shift @scores;
		$statesum -= shift @states;
	}
}
close IN;


monotonize();

my $index = 0;
$type1 = "";
$type2 = "";

open IN, "gzip -cd ".$root."_cooccur_pairs.tsv.gz |";

my $viruses = "viruses.species.tsv";
my %virus;
open VIR, "$viruses" or die "can't open $viruses for reading";
while (<VIR>) {
    chomp;
    my ($taxid, undef) = split /\t/;
    $virus{$taxid} = 1;
}

my %output;


while (<IN>) {
	s/\r?\n//;
	my ($newtype1, $identifier1, $newtype2, $identifier2, $score) = split /\t/;
	my $faketype1 = $newtype1;
	$faketype1 = 1 if $faketype1 < -1;
	my $faketype2 = $newtype2;
	$faketype2 = 1 if $faketype2 < -1;
	if ($faketype1 < $faketype2) {
		my $faketype3 = $faketype1;
		$faketype1 = $faketype2;
		$faketype2 = $faketype3;
	}
	$faketype2 = $faketype1 if $faketype2 > 0 and $faketype2 != $faketype1;
    $faketype1 = 9606;
    $faketype2 = 9606;
	if ($newtype1 ne $type1 or $newtype2 ne $type2) {
		$index = 0;
		$type1 = $newtype1;
		$type2 = $newtype2;
		if (exists $type_type_score_prob{$faketype1} and exists $type_type_score_prob{$faketype1}{$faketype2}) {
			@scores = sort {$b <=> $a} keys %{$type_type_score_prob{$faketype1}{$faketype2}};
		}
		else {
			@scores = ();
		}
	}
	my $prob = 0;
	next unless scalar(@scores);
	if ($score >= $scores[$index]) {
		$prob = $type_type_score_prob{$faketype1}{$faketype2}{$scores[$index]};
	}
	else {
		while ($index < $#scores and $score < $scores[$index+1]) {
			$index++;
		}
		if ($index == $#scores) {
			$prob = $type_type_score_prob{$faketype1}{$faketype2}{$scores[$index]};
		}
		else {
			my $score1 = $scores[$index];
			my $score2 = $scores[$index+1];
			my $prob1 = $type_type_score_prob{$faketype1}{$faketype2}{$score1};
			my $prob2 = $type_type_score_prob{$faketype1}{$faketype2}{$score2};
			$prob = ($prob1*($score-$score2)+$prob2*($score1-$score))/($score1-$score2);
		}
	}

    my $virus;
    next if exists $virus{$type1} and exists $virus{$type2} and $type1 != $type2; # skip v1-v2
    $virus = $type1 if exists $virus{$type1};
    $virus = $type2 if exists $virus{$type2};
    next if not $virus;


    my $print_type = $virus;
    if ($type1 != $type2) {
        $print_type = $type1 . "-" . $type2;
    }

    my $line = sprintf("cooccur\t%s\t%s\t%s\t%.3f\n", $print_type, $identifier1, $identifier2, $prob);
    push @{$output{$virus}}, $line;
    #printf OUT "cooccur\t%s\t%s\t%s\t%.3f\n", $print_type, $identifier1, $identifier2, $prob;
    
    if ($type1 != $type2) {
        $print_type = $type2 . "-" . $type1;
    }
    #printf OUT "cooccur\t%s\t%s\t%s\t%.3f\n", $print_type, $identifier2, $identifier1, $prob;
    $line = sprintf("cooccur\t%s\t%s\t%s\t%.3f\n", $print_type, $identifier2, $identifier1, $prob);
    push @{$output{$virus}}, $line;
}

my $c = scalar keys %output;

my $outfile = "string_virus_cooccur_interact";
open OUT, "| gzip -9 > $outfile.tsv.gz";
foreach my $taxid (keys %output) {
    foreach my $line (@{$output{$taxid}}) {
        print OUT $line;
    }
}


close IN;
close OUT;

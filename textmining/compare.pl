#!/usr/bin/perl -w

# compare db_... to virus.protein... output in different formats to ensure 


use strict;
use warnings;

my $file1 = shift;
my $file2 = shift;

open (my $fh1, '<', $file1) or die "can't open $file1 $!";
open (my $fh2, '<', $file2) or die "can't open $file2 $!";

my $taxid1 = "";
my $taxid2 = "";

my $line1 = <$fh1>; 
$line1 = <$fh1>; 
my $line2 = <$fh2>;
$line2 = <$fh2>;

my $acc = 0;
while (1) {
    $acc++;
    my @c1 = split /\t/, $line1;
    $taxid1 = $c1[1];
    my @c2 = split /\t/, $line2;
    my ($t, $p) = split /\./, $c2[1];
    $taxid2 = $t;

    if ($taxid1 != $taxid2) {
    print "acc $acc\n";
        print $c2[0] . "\t" . $c2[1] . "\n";
        $line2 = <$fh2>;
    }
    else {
        $line1 = <$fh1>;
        $line2 = <$fh2>;
    }
    last if not $line1;
}

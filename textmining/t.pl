use warnings;
use strict;

my $dir = {};
while (<>) {
    chomp;
    my @c = split "\t", $_;
    $dir->{$c[0]}->{$c[1]} += 1;
}

foreach my $key1 (keys %{$dir}) {
    foreach my $key2 (keys %{$dir->{$key1}}) {
        if (not defined $dir->{$key2} and not defined $dir->{$key2}->{$key1}) {
            print "no friend $key1 $key2\n";
        }
        if ($key1 eq $key2) {
            print "delete from links where entity1='$key1' and entity2='$key2'\n";
        }
    }
}

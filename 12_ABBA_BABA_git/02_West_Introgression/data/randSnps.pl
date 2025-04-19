#!/usr/bin/perl
use strict;
use warnings;
my %snps_per_seq = ();
while (<>) {
    if (/^#/) { # print headers as the original
        print; 
    }
    else {
        my ($seq_id) = split (/\t/, $_);
        push @{ $snps_per_seq{ $seq_id } }, $_; # store a hash per sequence, with possitions as array
    }
}
foreach my $seq_id (sort keys %snps_per_seq) {
    my @snps = @{ $snps_per_seq{ $seq_id } };
    print $snps[ int rand @snps ]; # grab one random position from the array
}

#usage: perl randSnps.pl < input.vcf > output.vcf



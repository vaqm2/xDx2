#!/usr/bin/env perl

use strict;
use warnings;
use IO::File;

my $fh = IO::File->new("$ARGV[0]") || die "Cannot open file: $ARGV[0]!\n";
my $map = {};

while(my $line = $fh->getline) {
    chomp($line);
    my @lineContents = split(/\t/, $line);
    my $raw_id       = $lineContents[0];
    my $hgnc_symbol  = $lineContents[5];
    $map->{$raw_id}  = $hgnc_symbol;
}

$fh->close;

$fh = IO::File->new("$ARGV[1]") || die "Cannot open file: $ARGV[1]!\n";

while(my $line = $fh->getline) {
    chomp($line);
    if($line =~ /^\#/) {
        print $line."\n";
        next;
    }
    else {
        my @lineContents = split(/\s+/, $line);
        my $raw_id       = $lineContents[0];

        if(exists $map->{$raw_id}) {
            print $map->{$raw_id};

        for my $index(1..$#lineContents) {
            print " ".$lineContents[$index];
        }

        print "\n";
        }
        else {
            next;
        }
    }
}

$fh->close;
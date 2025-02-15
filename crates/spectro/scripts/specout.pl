#!/usr/bin/perl

use v5.36;

# translate spectro output geometries to normal atoms

# default masses to atom labels
my %atoms = (
    "1.0078250" => "H",
    "12.0000000" => "C",
    "14.0030740" => "N",
    "18.9984032" => "F",
    "23.9850423" => "Mg",
);

while (<>) {
    my @split = split " ";
    $split[4] = $atoms{$split[4]};
    printf "%2s%18.10f%18.10f%18.10f\n", $split[4], $split[1], $split[2], $split[3];
}

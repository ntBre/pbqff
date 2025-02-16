#!/usr/bin/perl

use v5.36;

# translate spectro input geometries to normal atoms

my %atoms = (
    "1.00" => "H",
    "2.00" => "He",
    "3.00" => "Li",
    "4.00" => "Be",
    "5.00" => "B",
    "6.00" => "C",
    "7.00" => "N",
    "8.00" => "O",
    "9.00" => "F",
);

while (<>) {
    my @split = split " ";
    $split[0] = $atoms{$split[0]};
    printf "%2s%18.10f%18.10f%18.10f\n", $split[0], $split[1], $split[2], $split[3];
}

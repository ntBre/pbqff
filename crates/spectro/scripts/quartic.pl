#!/usr/bin/perl

use v5.36;

# maybe a script for generating annoying utility implementations for Quartic

my @display;
my @subtract;
while(<>) {
    my $display = s/\s+pub\(crate\) (.*): f64,/$1/r;
    my $subtract = s/\s+pub\(crate\) (.*): f64,/$1: self.$1 - rhs.$1,/r;
    push @subtract, $subtract;
}

say "Self {";
for (@subtract) {
    print;
}
say "}";

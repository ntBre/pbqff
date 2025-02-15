#!/usr/bin/perl

use v5.36;

my @lines;
while(<>) {
    push @lines, $_;
}

my $text = join "", @lines;

$text =~ s/&nalgebra::Matrix<
\s*f64,
\s*nalgebra::Dynamic,
\s*nalgebra::Const<1>,
\s*nalgebra::VecStorage<f64, nalgebra::Dynamic, nalgebra::Const<1>>,
\s*>,
/\&Dvec,\n/g;

$text =~ s/nalgebra::Matrix<
\s*f64,
\s*nalgebra::Dynamic,
\s*nalgebra::Dynamic,
\s*nalgebra::VecStorage<f64, nalgebra::Dynamic, nalgebra::Dynamic>,
\s*>,
/Dmat,\n/g;

say "$text";

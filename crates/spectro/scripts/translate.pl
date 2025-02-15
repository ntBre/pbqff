#!/usr/bin/perl

use warnings;
use strict;

sub work {
  my $hold = shift;
  $hold =~ m#(\d+)# ? $hold : "=$hold"
}

while (<>) {
  # add semi-colons
  s/$/;/;

  # lowercase letters
  s/(\w)/\L$1/g;

  # named numbers
  s/zero/0.0/g;
  s/one/1.0/g;
  s/two/2.0/g;

  # scientific notation
  s/\.0d0/.0/g;

  # simple ident followed by =
  s/^\s*([a-z][0-9])+\s*=(.*)$/let $1 = $2/;

  # do loop
  s/do\s+\d+\s+(\w)=(\w),(\w).*/
    sprintf "for $1 in %d..%s {", $2-1, work($3)/ex;

  # continue
  s/.*continue.*/}/g;

  # ident followed by matrix subscripts
  s/([a-z0-9]+)(\([^)]+\))/$1\[$2\]/g;

  # square
  s/\*\*2/.powi(2)/g;

  # abs and trig
  s/dabs/f64::abs/g;
  s/dcos/f64::cos/g;
  s/dsin/f64::sin/g;

  s/abs\[\(([^)]+)\)\]/abs($1)/g;
  s/cos\[\(([^)]+)\)\]/cos($1)/g;
  s/sin\[\(([^)]+)\)\]/sin($1)/g;

  print;
}

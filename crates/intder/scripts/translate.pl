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

  # simple ident followed by =
  s/^\s*([a-z][0-9])+\s*=(.*)$/let $1 = $2/;

  # do loop
  s/do\s+\d+\s+(\w)=(\w),(\w).*/
    sprintf "for $1 in %d..%s {", $2-1, work($3)/ex;

  # continue
  s/.*continue.*/}/g;

  # ident followed by matrix subscripts
  s/([a-z0-9]+)(\([^)]+\))/$1\[$2\]/g;

  print;
}

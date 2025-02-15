#!/usr/bin/awk -f

# for extracting matrices from intder.out. currently hard-coded for A MATRIX

/A MATRIX/ { print "\n"; next }
$1 ~ /I=/ { print ""; next }
{
    for (i = 1; i <= NF; i++) printf "%10.6f", $i
}

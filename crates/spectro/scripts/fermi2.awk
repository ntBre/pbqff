#!/usr/bin/awk -f

{
    printf "Fermi2::new(%5d,%5d,%5d),\n", $2-1, $3-1, $4-1
}

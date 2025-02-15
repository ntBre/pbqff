#!/usr/bin/awk -f

{
    for (i = 1; i < NF; i++) {
	printf "%s,", $i
    }
    printf "%s;\n", $NF
}

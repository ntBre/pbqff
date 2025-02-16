#!/usr/bin/awk -f

# print an integer matrix from Fortan, subtracting 1 from each element
{
    for (i = 1; i < NF; i++) {
	printf "%5d,", $i-1
    }
    printf "%5d;\n", $NF-1
}

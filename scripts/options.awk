#!/usr/bin/awk -f

/struct Args/,/^}$/ {
    if ($1 ~ /\/\/\//)  {
	for (i = 2; i <= NF; i++) {
	    buf = buf " " $i
	}
    } else if (/^$/ || /struct/ || /^}$/ || $1 ~ /^#/) {
	# pass
    } else { # field name
	field = substr($1, 1, length($1) - 1) # trim colon
	type = substr($2, 1, length($2) - 1) # trim comma
	if (field !~ /infile|json/)
	    printf ".TP\n.B %s \\fI%s\\fR\n%s\n", field, type, substr(buf, 2)
	buf = ""
    }
}

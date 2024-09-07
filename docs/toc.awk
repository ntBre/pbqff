#!/usr/bin/awk -f

BEGIN { print "<ul>" }
END { print "</ul>" }

/head[23]\(/ {
	printf "\t"
}

/head2\(/ {
	if (sublist) {
		print "</ul>" # close previous sublist
		sublist = 0
	}
	match($0, /head2\(([^)]+)\)/, arr)
	printf "<li><a href=\"#%s\">%s</a></li>\n", make_id(arr[1]), arr[1]
}

/head3\(/ {
	if (!sublist) {
		print "<ul>"
	}
	sublist = 1
	printf "\t"
	match($0, /head3\(([^)]+)\)/, arr)
	printf "<li><a href=\"#%s\">%s</a></li>\n", make_id(arr[1]), arr[1]
}

function make_id(s,   ret) {
	ret = tolower(s)
	gsub(/ +/, "-", ret)
	return ret
}

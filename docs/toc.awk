#!/usr/bin/awk -f

BEGIN { print "<ul>" }
END { print "</ul>" }

/head[234]\(/ {
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

/head4\(/ {
	print "<ul>" # always add a list level for this, assume it's alone
	printf "\t\t"
	match($0, /head4\(([^)]+)\)/, arr)
	printf "<li><a href=\"#%s\">%s</a></li>\n", make_id(arr[1]), arr[1]
	print "</ul>"
}

function make_id(s,   ret) {
	ret = tolower(s)
	gsub(/ +/, "-", ret)
	return ret
}

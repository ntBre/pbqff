#!/usr/bin/awk -f

function make_id(s,   ret) {
	ret = tolower(s)
	gsub(/ /, "-", ret)
	return ret
}

BEGIN {
	HEAD = "head([2345])\\(([^)]+)\\)"
}
END { print "</ul>" }

$0 ~ HEAD {
	match($0, HEAD, arr)
	cur = arr[1] # current heading level
	while (p > 0 && cur < buf[p-1]) {
		#printf("closing list with cur = %d, last = %d\n", cur, buf[p-1])
		print "</ul>"
		p--
	}
	if (buf[p-1] < cur) {
		#printf "opening new list with cur = %d, last = %d\n", cur, buf[p-1]
		print "<ul>"
	}
	if (buf[p-1] != arr[1])
		buf[p++] = arr[1] # heading level stack

	printf "<li><a href=\"#%s\">%s</a></li>\n", make_id(arr[2]), arr[2]
}

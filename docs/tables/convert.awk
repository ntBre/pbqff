#!/usr/bin/env -S awk --csv -f

# generate the body of an HTML table from a CSV file

NR == 1 {
	print "<tr>"
	for (i = 1; i <= NF; i++)
		printf "<th>%s</th>\n", $i
	print "</tr>"
}

NR > 1 {
	print "<tr>"
	for (i = 1; i <= NF; i++)
		printf "<td>%s</td>\n", $i
	print "</tr>"
}

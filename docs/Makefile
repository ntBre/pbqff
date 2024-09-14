toc.html: template.html toc.awk
	./toc2.awk $< > $@

index.html: template.html toc.html
	m4 $< > $@

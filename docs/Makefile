toc.html: template.html toc.awk
	./toc.awk $< > $@

index.html: template.html toc.html
	m4 $< > $@

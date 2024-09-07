toc.html: template.html
	./toc.awk $^ > $@

index.html: template.html toc.html
	m4 $< > $@

%.html : %.md
	pandoc $< > $@

%.html : %.Rmd
	R --vanilla -e "require(rmarkdown);render(\"$<\",output_file=\"$@\")"

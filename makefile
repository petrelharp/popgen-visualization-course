%.html : %.md
	pandoc $< > $@

%.html : %.Rmd
	R --vanilla -e "require(rmarkdown);render(\"$<\",output_file=\"$@\")"

slides.html : slides.Rmd gideons-ibde-plot-1.png gideons-ibde-plot-2.png
	R --vanilla -e "require(rmarkdown);render(\"$<\",output_file=\"$@\")"

gideons-ibde-plot-%.pdf : gideons-ibde-plot.ink.svg
		./scripts/export-layers-svg.sh $< ibd0 ibd$* >$@

%.png : %.pdf
	convert -density 300 $< -flatten $@


ifeq ($V,1)
silent=
silent_stdout=
else
silent=@
silent_stdout= > /dev/null
endif

docdvi: tfm.tex
	cp tfm.tex tfm.bak
	touch --reference=tfm.tex tfm.bak
	(printf "%s" "\def\fixedpdfdate{"; date +'D:%Y%m%d%H%M%S%:z' -d @$$(stat --format=%Y tfm.tex) | sed "s/:\([0-9][0-9]\)$$/'\1'}/g") > tfm-deterministic.tex
	printf "%s\n" "\pdfinfo{" >> tfm-deterministic.tex
	printf "%s\n" "  /CreationDate (\fixedpdfdate)" >> tfm-deterministic.tex
	printf "%s\n}\n" "  /ModDate (\fixedpdfdate)" >> tfm-deterministic.tex
	cat tfm.tex >> tfm-deterministic.tex
	mv tfm-deterministic.tex tfm.tex
	touch --reference=tfm.bak tfm.tex
	touch tfm.ind
	latex tfm $(silent_stdout)
	latex tfm $(silent_stdout)
	makeindex tfm
	latex tfm $(silent_stdout)

tfm.pdf: docdvi
	latex tfm $(silent_stdout)
	pdflatex tfm $(silent_stdout)
	sed -b -i 's,^/ID \[.*\]$$,/ID [<0> <0>],g' tfm.pdf
	mv tfm.bak tfm.tex

docs: tfm.pdf
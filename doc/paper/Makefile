PAPER = main.pdf
TEXFILES = *.tex
BIBFILES = references.bib

all: pdf

clean:
	-rm -f *.aux *.log *.dvi *.blg *.bbl *.inx *.aut *.toc *.out *Notes.bib $(PAPER)

pdf: $(PAPER)

main.pdf: $(BIBFILES) $(TEXFILES)
	pdflatex main.tex
	pdflatex main.tex
	bibtex   main
	pdflatex main.tex
	pdflatex main.tex

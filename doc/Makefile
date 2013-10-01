                                                                                                                                                               
all: pdf

clean:
	-rm -f *.aux *.log *.dvi *.blg *.out ms1.ps

# ----------------------------------------------------------------------------

pdf: ms1.pdf

ps: ms1.ps

# ----------------------------------------------------------------------------

# Use scriptutils "bubble" script to generate .bbl from .bib database:
# ms.bbl: ms.tex references.bib
# 	bubble -f ms.tex references.bib >& ms.bbl.log 
#       
# ms.dvi: ms.bbl ms.tex
# 	latex ms.tex
# 	latex ms.tex
      
ms1.dvi: master.bib ms1.tex
	latex ms1.tex
	latex ms1.tex
# 	bibtex ms1
# 	latex ms1.tex
# 	latex ms1.tex

ms1.pdf: ms1.dvi
	\dvipdf -smsSIZE=letter ms1.dvi

ms1.ps: ms1.dvi
	\dvips -t letter -o ms1.ps ms1.dvi

# ----------------------------------------------------------------------------

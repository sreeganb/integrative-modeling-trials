# Specify the LaTeX document and the PDF output file
DOCUMENT = tips.tex
PDF_FILE = tips.pdf
BIB_FILE = tips.bib

# Rule for compiling the LaTeX document
$(PDF_FILE): $(DOCUMENT) $(BIB_FILE)
	pdflatex $(DOCUMENT)
	bibtex $(basename $(DOCUMENT))
	pdflatex $(DOCUMENT)
	pdflatex $(DOCUMENT)

# Rule for cleaning up temporary files
clean:
	rm -f *.aux *.log *.dvi *.lof *.lot *.toc *.blg *.bbl *.pdf
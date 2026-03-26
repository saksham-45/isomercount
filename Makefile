.PHONY: ass1 ass2 pdfs

ass1:
	cd docs && pdflatex -interaction=nonstopmode -halt-on-error -jobname ass1 cyclic_counting_proposal.tex
	cd docs && bibtex ass1
	cd docs && pdflatex -interaction=nonstopmode -halt-on-error -jobname ass1 cyclic_counting_proposal.tex
	cd docs && pdflatex -interaction=nonstopmode -halt-on-error -jobname ass1 cyclic_counting_proposal.tex
	mkdir -p output/pdf
	cp -f docs/ass1.pdf output/pdf/ass1.pdf

ass2:
	cd docs && pdflatex -interaction=nonstopmode -halt-on-error -jobname ass2 snake_polyomino_proposal.tex
	cd docs && bibtex ass2
	cd docs && pdflatex -interaction=nonstopmode -halt-on-error -jobname ass2 snake_polyomino_proposal.tex
	cd docs && pdflatex -interaction=nonstopmode -halt-on-error -jobname ass2 snake_polyomino_proposal.tex
	mkdir -p output/pdf
	cp -f docs/ass2.pdf output/pdf/ass2.pdf

pdfs: ass1 ass2

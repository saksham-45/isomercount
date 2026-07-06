.PHONY: ass1 ass2 pdfs final_report cyclic_report snake_report reports

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

final_report:
	cd docs/final_report && pdflatex -interaction=nonstopmode -halt-on-error report.tex
	cd docs/final_report && pdflatex -interaction=nonstopmode -halt-on-error report.tex
	mkdir -p output/pdf
	cp -f docs/final_report/report.pdf output/pdf/final_report.pdf

cyclic_report:
	cd docs/final_report && pdflatex -interaction=nonstopmode -halt-on-error cyclic_report.tex
	cd docs/final_report && pdflatex -interaction=nonstopmode -halt-on-error cyclic_report.tex
	mkdir -p output/pdf
	cp -f docs/final_report/cyclic_report.pdf output/pdf/cyclic_report.pdf

snake_report:
	cd docs/final_report && pdflatex -interaction=nonstopmode -halt-on-error snake_report.tex
	cd docs/final_report && pdflatex -interaction=nonstopmode -halt-on-error snake_report.tex
	mkdir -p output/pdf
	cp -f docs/final_report/snake_report.pdf output/pdf/snake_report.pdf

reports: cyclic_report snake_report

pdfs: ass1 ass2 final_report cyclic_report snake_report

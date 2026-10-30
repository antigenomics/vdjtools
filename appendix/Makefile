# Compile the recombination-model appendix. Requires lualatex (OldStandard OpenType fonts;
# pinned in .latexmkrc) + bibtex and the two frozen network PDFs.
murugan_model.pdf: murugan_model.tex refs.bib bn_trb.pdf bn_trd_dd.pdf .latexmkrc
	latexmk -lualatex -interaction=nonstopmode -halt-on-error murugan_model.tex

.PHONY: clean
clean:
	latexmk -C murugan_model.tex

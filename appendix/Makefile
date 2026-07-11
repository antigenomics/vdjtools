# Compile the recombination-model appendix. Requires lualatex (OldStandard OpenType fonts;
# pinned in .latexmkrc) + bibtex. No external assets.
murugan_model.pdf: murugan_model.tex refs.bib
	latexmk -lualatex -interaction=nonstopmode -halt-on-error murugan_model.tex

.PHONY: clean
clean:
	latexmk -C murugan_model.tex

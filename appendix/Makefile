# Compile the recombination-model appendix. Requires lualatex (OldStandard OpenType fonts;
# pinned in .latexmkrc) + bibtex and the two frozen network PDFs.
vdj_rearrangement_model.pdf: vdj_rearrangement_model.tex refs.bib bn_trb.pdf bn_trd_dd.pdf .latexmkrc
	latexmk -lualatex -interaction=nonstopmode -halt-on-error vdj_rearrangement_model.tex

.PHONY: clean
clean:
	latexmk -C vdj_rearrangement_model.tex

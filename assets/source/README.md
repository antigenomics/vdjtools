# README analysis atlas

Build with `python3 assets/source/build_analysis_atlas.py`, then render with
`rsvg-convert assets/analysis-atlas-light.svg -o assets/analysis-atlas.png`.
The palette follows isalgo.dev; diagrams other than the motif panel are schematic.

## Recorded motif example

Source: [VDJdb 2025-02-21 release](https://github.com/antigenomics/vdjdb-db/releases/tag/2025-02-21),
`motif_pwms.txt` and `cluster_members.txt`, cluster `H.B.NLVPMVATV.10`,
human TRB, NLVPMVATV, 13-residue anchor-inclusive junctions, 9 member sequences.
These are public released observations; no model fitting or analysis was rerun.

`motif-example.tsv` retains the released `pos`, `aa`, `freq` and `freq.bg`
columns verbatim. All 13 motif frequency columns sum to 1 (checked).
The logo displays signed `freq - freq.bg` for the observed residues,
with positive letters above the baseline and negative letters below it.
The background is conditioned on V, J and junction length; it discounts
conserved germline contributions. Background-only residues absent from the
released motif rows are not drawn. This is a frequency-difference logo,
not an information-content or log-odds logo.

`kmer-example.tsv` counts overlapping 3-mers after dropping four residues
from each end of each of those nine junctions. The cloud displays the seven
most frequent words; font size follows their counts. Fixed four-residue trimming
is a display convention, not an inferred germline boundary.

To recover the original tables, download the release archive:

```sh
curl -L https://github.com/antigenomics/vdjdb-db/releases/download/2025-02-21/vdjdb-2025-02-21.zip -o /tmp/vdjdb-2025-02-21.zip
```

VDJdb data attribution: VDJdb Developers and Maintainers; the release includes
its CC BY-ND 4.0 license. The TSV files retain the selected source frequencies.

The amino-acid strip is the per-position maximum-frequency consensus of the
recorded motif PWM (not necessarily one observed junction). Its tile colors use
`hydropathy` from the package's `aa_property_table.txt`, linearly mapped from
−4.5 (hydrophilic, teal) to +4.5 (hydrophobic, amber).

PWM letter colors encode residue classes independently of enrichment sign:
G amber, P brown, C ochre, ST teal, NQ violet, DE red, KRH blue,
FWY mauve, and AILMV slate. Light/dark variants preserve these assignments;
letter position above/below the baseline continues to encode the signed difference.

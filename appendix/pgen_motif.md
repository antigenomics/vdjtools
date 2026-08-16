# Memo: Pgen of a motif — expose the masked DP, then (optionally) an offset automaton

2026-08-16. Requested for the mhcmatch immunogenicity / precursor-frequency work
(`~/vcs/projects/2026-mhcmatch-benchmark`). Status: **Proposal 1 implemented**
(`native.pgen_aa_degenerate` / `pgen_aa_degenerate_batch`, `tests/python/test_native_degenerate.py`);
Proposal 2 (the offset automaton) is still a proposal.

## The ask

Given an epitope, we derive a **paratope motif** — either a VDJdb cluster PWM
(`motif_pwms.txt`, V/J/length-pinned) or a set of V-pinned k-mers (V+3-mer, V+gapped-4-mer) — and we
need its **generation probability**: the total Pgen of every junction the motif matches. That number
is the model-side estimate of an epitope's T-cell precursor frequency.

## Finding: the engine already exists and is already validated

`pgen_aa_masked(const PackedModel&, const uint64_t* allowed, int L, int v_idx, int j_idx)`
— `src/pgen.cpp:701` — is exactly a **degenerate-sequence DP**. `allowed[c]` is a 64-bit codon
bitmask per position, so a position can permit any residue subset at no extra cost: the transfer
matrix contracts over the allowed codons either way.

It is not new or untested code. It is the engine behind both public entry points:

- `pgen_aa` (`src/pgen.cpp:711`) — every mask is a single residue (`mask_for_aa`).
- `pgen_aa_hamming1` (`src/pgen.cpp:721`) — sets one position to `mask_wildcard()` at a time and
  applies inclusion–exclusion. Pinned against brute-force enumeration of all 19L neighbours at
  rel 1e-9 (`tests/python/test_overlap_alice.py:33`) and against OLGA at rtol 1e-9
  (`tests/python/test_native.py:91`).

So the arbitrary-subset case is exercised on every Hamming-1 call today; only the *wildcard* subset
is reachable from Python.

**The gap is purely a binding.** `pgen_aa_masked` sits in an anonymous namespace
(`src/pgen.cpp:700`, `}  // namespace` at :709), and `src/_bindings.cpp` exposes only `pgen_nt`,
`pgen_aa`, `pgen_aa_hamming1`, `pgen_aa_batch` (lines 64/68/72/76). Nothing else is missing.

Correction from the implementation: promoting it is a **two**-place edit, not one. `pgen_nt`'s
in-frame fast path calls `pgen_aa_masked` through a *second*, forward declaration in the file's
first anonymous namespace (`src/pgen.cpp:188`). Leave that in place and the promoted symbol makes
every call site ambiguous — the build fails at six of them. Delete the forward declaration; the
header declaration covers those calls.

## Proposal 1 — expose the masked DP (small, unlocks most of the use case)

1. Promote `pgen_aa_masked` out of the anonymous namespace; declare it in
   `include/vdjtools/model.hpp` beside `pgen_aa_hamming1` (:65).
2. Bind it in `src/_bindings.cpp` taking a per-position list of allowed residues rather than raw
   `uint64_t` — build the mask C++-side by OR-ing `mask_for_aa` over each position's residue set, so
   the codon encoding never leaks into Python.
3. Python wrapper in `model/native.py` next to `pgen_aa`:

   ```python
   pgen_aa_degenerate(model, allowed, v=None, j=None) -> float
   # allowed: sequence of length L, each item a str of permitted residues
   #   "C" pins,  "ILVF" allows a subset,  "" or "X" = wildcard
   ```

   Plus a batch form mirroring `pgen_aa_batch(threads=0)`.

**Invariants worth asserting in the test**, all free consequences of the identity:

- `pgen_aa_degenerate(m, list(seq)) == pgen_aa(m, seq)` exactly.
- All-wildcard equals `P(L, V, J)`, the model's length marginal.
- Monotone under set inclusion: widening any position never decreases Pgen.
- The Hamming-1 identity reproduces `pgen_aa_hamming1` from the public API alone.

### What this unlocks immediately

- **VDJdb cluster PWMs.** `motif_pwms.txt` is V/J/length-pinned (`v.segm.repr`, `j.segm.repr`,
  `len`), so a motif is precisely a per-position residue set. Thresholding `I`/`freq` per position
  gives `allowed`; one masked call gives the cluster's Pgen. No automaton, no enumeration.
- **Any wildcard/degenerate query**, including the Hamming-1 ball at explicit positions.
- **Gapped patterns at a fixed offset** — a gap is a wildcard position.

Note `mask_for_aa` (`src/pgen.cpp:268`) matches by exact character against the genetic code, so `X`
currently yields an **empty** mask and scores 0.0, not a wildcard. The wrapper must map `X`/`""` to
`mask_wildcard()` explicitly, or `pgen_aa(m, "…X…") == 0.0` becomes a silent trap.

## Proposal 2 — "junction contains this k-mer at any offset" (larger, defer)

Proposal 1 handles anything pinned to fixed positions. It does **not** handle an unanchored k-mer,
where the same k-mer may occur at several offsets and lengths; summing over offsets double-counts
junctions matching twice.

The exact object is the complement:

```
Pgen(contains ≥1 of Π)  =  1 − Pgen(avoids every π ∈ Π)
```

Build the **Aho–Corasick automaton** over `Π`, delete accepting states, and run the recombination DP
over the **product** of the packed model with that automaton. Overlap, multiplicity and offset are
absorbed into automaton state — one pass, cost linear in `|model states| × Σ|π|`, no
inclusion–exclusion. To require ≥ j distinct matches, add a saturating counter to the state.

Since V is a choice node in the model, `Pgen(V=v, contains π)` is just the per-V DP — so **V-pinned
k-mers are cheaper than unpinned ones**, which is the usual feature spec anyway.

**Recommendation: do Proposal 1 first and only.** For V/J/length-pinned VDJdb cluster motifs the
masked DP is exact and sufficient; the automaton is worth building only if unanchored k-mer sets
turn out to be the better epitope→paratope representation. Decide that empirically, not now.

## Traps carried over

- **Junction, not CDR3.** Anchors included; an IMGT CDR3 scores 0.0 silently.
- **Allele-resolution V/J.** `native._gene_idx` raises on a gene-level name by design; with
  `collapse=True` (the default) pass `TRBV27*01`.
- **Marginal vs conditioned.** `pgen_aa(m, seq)` marginalises over V/J; `pgen_aa(m, seq, v, j)` is a
  joint. VDJdb motifs are V/J-pinned, so the conditioned form is the right one — do not mix them in
  one comparison.
- **Mouse:** `source="arda"` only, TRA/TRB only.

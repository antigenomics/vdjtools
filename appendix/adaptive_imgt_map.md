# Adaptive/immunoSEQ → IMGT gene names: a CDR-validated table

**Status: applied in 3.2.0 (2026-08-09).** `io/convert.py::_adaptive_to_imgt` is the §5 patch;
`tests/python/test_convert.py::test_adaptive_to_imgt` is the §7 self-check as a unit test. This note
is kept as the rationale and the regeneration recipe for the shipped table
`python/vdjtools/resources/adaptive_imgt_map.tsv` (builder: `appendix/build_adaptive_imgt_map.py`).
Diagnosis below is as of 2026-08-08, i.e. against the pre-fix implementation.

---

## 1. What is wrong today

`vdjtools/io/convert.py`:

```python
_ZERO_PAD = re.compile(r"0([1-9])")

def _adaptive_to_imgt(field):
    gene = extract_vdj(field)                       # first comma-tie, allele after "*" stripped
    if gene is None or gene.lower() == "unresolved":
        return None
    return _ZERO_PAD.sub(r"\1", gene.replace("TCR", "TR"))
```

The zero-strip is a **global** substitution, so the trailing `-01` of an Adaptive token is always
re-emitted as an IMGT **subgroup**. For most tokens that trailing group is an **allele**, and the
result is a gene name that does not exist.

Measured on the 161 distinct Adaptive tokens of the IMMREP25 release + the pairSEQ mock cohort
(44,000 gene calls):

| | tokens | gene calls |
|---|---|---|
| `_adaptive_to_imgt` agrees with the validated map | 61 | 21,942 |
| `_adaptive_to_imgt` disagrees | **100** | **22,058** |
| `_adaptive_to_imgt` output absent from the IMGT human reference | **100** | **22,058** |

Every single disagreement produces a **non-existent gene name** — i.e. every one of these is a
straight defect, not a nomenclature preference. Four failure modes:

**C. Allele suffix read as a subgroup** — 86 tokens, 19,579 calls. The large one.

| Adaptive | `_adaptive_to_imgt` | correct | note |
|---|---|---|---|
| `TCRAJ39-01` | `TRAJ39-1` | `TRAJ39` | **no** human TRAJ gene has a subgroup (61 TRAJ genes, zero `-N`) |
| `TCRAV22-01` | `TRAV22-1` | `TRAV22` | TRAV22 is subgroup-less |
| `TCRBV09-01` | `TRBV9-1` | `TRBV9` | 795 calls |
| `TCRBV27-01` | `TRBV27-1` | `TRBV27` | 671 calls |
| `TCRBV02-01` | `TRBV2-1` | `TRBV2` | 482 calls |
| `TCRBD01-01` | `TRBD1-1` | `TRBD1` | D genes too — `_adaptive_call` feeds `D_CALL`; not among the 161 (no paired cohort carries D calls), added to the table via `--also` |

…and **the same textual pattern is correct** for `TCRAV01-01 → TRAV1-1`, `TCRBV29-01 → TRBV29-1`,
`TCRBV25-01 → TRBV25-1`, `TCRAV08-03 → TRAV8-3`.

**A. Adaptive ties are not split** — 3 tokens, 410 calls. `extract_vdj` splits on `,` only, so a
slash tie survives into the output: `TCRBV03-01/03-02 → TRBV3-1/3-2`, `TCRBV06-02/06-03 →
TRBV6-2/6-3`, `TCRBV12-03/12-04 → TRBV12-3/12-4`.

**B. Family-only calls pass through verbatim** — 6 tokens, 1,090 calls.
`TCRBV20-X → TRBV20-X` (680 calls), `TCRBV12-X → TRBV12-X`, `TCRBV05-X`, `TCRBV06-X`, `TCRBV07-X`,
`TCRBV11-X`.

**D. Co-locus (`/DVn`) names are not reconstructed** — 5 tokens, 979 calls.
`TCRAV38-02 → TRAV38-2`, but IMGT calls it `TRAV38-2/DV8`. Likewise `TCRAV14-01 → TRAV14/DV4`,
`TCRAV23-01 → TRAV23/DV6`, `TCRAV29-01 → TRAV29/DV5`, `TCRAV36-01 → TRAV36/DV7`.

### Consequence in practice

Any consumer that resolves gene names against a germline reference silently loses the rows. The
concrete case that surfaced this: tcrdist3 filters `v_a_gene`/`j_a_gene`/`v_b_gene`/`j_b_gene`
against its own germline set and dropped **100 %** of both affected cohorts:
`input rows=1000; kept=0 (dropped 1000 unknown-gene)` for the IMMREP25 positives and the same for
the pairSEQ mock. With the validated map: `kept=1000 (dropped 0)` and `kept=993 (dropped 7)` — the 7
are `TCRDV01-01` = `TRDV1`, a delta-locus V that rearranges into the alpha chain and that tcrdist3's
alphabeta frame will not accept in the `v_a` slot (it files `TRDV1` under chain `B`). The mapping is
correct; the reference's slot rule is what rejects it.

## 2. Why a regex cannot get this right

Whether the trailing group of `TCRxYNN-MM` is a **subgroup** or an **allele** is a property of the
gene family, not of the string:

* `TRAV1` has subgroups (`TRAV1-1`, `TRAV1-2`) → `TCRAV01-01` = subgroup 1.
* `TRAV22` has none → `TCRAV22-01` = allele 01 of `TRAV22`.
* No human `TRAJ` has a subgroup (61 genes, 68 alleles, zero `-N`), so `-01` is *always* an allele
  there — 54 of the 161 tokens, one rule.
* `TRBV9`, `TRBV27`, `TRBV28`, `TRBV19`, `TRBV2`, `TRBV13`–`TRBV18`, `TRBV30`, `TRBD1`, `TRBD2` have
  none; `TRBV29-1`, `TRBV25-1`, `TRBV20-1`, `TRBV3-1`, `TRBV5-1` do.

The two readings are textually identical. The decision needs the **IMGT gene list**, and the choice
among alleles/subgroups that the gene list leaves open needs **sequence** evidence. So the fix is a
table, not a better pattern.

## 3. How the table was built (the CDR-matching algorithm)

`appendix/build_adaptive_imgt_map.py`, stdlib only. Per Adaptive token:

1. **Parse** `TCR<locus><segment><family>[-<group>]`, `<group>` ∈ `NN` | `X` | `orNN_MM`. Split
   slash ties into one call each.
2. **Enumerate candidates** from the IMGT gene list, offering *both* readings — the subgroup reading
   `TR{loc}{seg}{fam}-{group}` **and** the allele reading `TR{loc}{seg}{fam}` — and keeping whichever
   exist. `/DVn` and `/ORn-n` are stripped to get a parseable core, so `TRAV14/DV4` is reachable from
   `TCRAV14-01`; orphons (`/OR`) are excluded unless the token itself says `orNN_MM`. `-X` expands to
   every gene of the family.
3. **Select by germline CDR sequence.** The IMMREP25 release ships, for every record, the germline
   `tcra_cdr1`/`tcra_cdr2` (and the `tcrb_*` twins) of the V gene it called. Ungap tcrdist3's
   germline `cdrs` field (`CDR1;CDR2;CDR2.5;CDR3-start`) and keep the candidate alleles whose CDR1
   **and** CDR2 equal the pair the release itself reports for that token. This validates the mapping
   against the source data instead of asserting it. Boundaries agree exactly — e.g. `TCRAV01-02`
   ships `TSGFNG`/`NVLDGL`, and `TRAV1-2*01` is `TSG......FNG`/`NVL....DGL`.
4. **J genes have no CDR1/CDR2.** The Adaptive J family number is unique in IMGT (at most one gene
   per family number), so the gene list is the primary evidence and the CDR3 3′ end is the
   cross-check. Only the germline **3′** end of a J survives recombination — the 5′ end is nibbled
   and overwritten by junctional insertions — so score candidates by mean **common-suffix length**
   between the observed CDR3s and the germline J segment, and accept at ≥ 5 aa (the Phe/Trp118
   anchor plus four). A whole-segment match is the wrong test and rejects almost every real J.
5. **Ambiguity is recorded, not hidden.** Alleles (or genes) sharing an identical CDR1+CDR2 are
   indistinguishable from this evidence; the lowest-numbered allele is chosen and
   `allele_ambiguous=True` is set. Because tcrdist3 and any CDR-based metric only ever read
   CDR1/CDR2/CDR2.5 from the V gene, a CDR-identical substitute is a **lossless** choice there.
6. **Nothing is invented.** If no candidate reproduces the observed CDRs, the token is left
   unmapped with `status=unresolved_dropped`.

### Result

163 tokens, none unmapped — the 161 observed in the two cohorts plus `TCRBD01-01`/`TCRBD02-01`,
added via `--also` because `read_immunoseq` maps `D_CALL` too but no paired cohort carries a D call:

| status | tokens | records |
|---|---|---|
| `validated_by_cdr` | **138** | 39,285 |
| `gene_list_unique` | 16 | 3,215 |
| `tie_resolved` | 9 | 1,500 |
| `unresolved_dropped` | 0 | 0 |
| `no_candidate` | 0 | 0 |

`gene_list_unique` = the gene list is unambiguous but the inputs carry no usable sequence evidence:
2 `--also` D tokens (`TCRBD01-01 → TRBD1`, `TCRBD02-01 → TRBD2`); 6 V tokens seen only in the
pairSEQ mock, which ships no CDR columns (`TCRAV14-01`, `TCRAV41-01`, `TCRDV01-01`, `TCRBV03-01`,
`TCRBV12-03`, `TCRBV23-01`); 6 J tokens with too few records to score (`TCRAJ01-01` n=2,
`TCRAJ02-01` n=1, `TCRAJ14-01` n=3, `TCRAJ19-01` n=1, `TCRAJ46-01` n=1, `TCRAJ60-01` n=1); and
`TCRBJ01-02` / `TCRBJ02-07`, whose 3′ agreement (4.90 and 4.83 aa over 621 and 2,496 rows) sits just
under the 5 aa bar because their germline segments are only 6 aa long — both are certainly right,
the bar is simply not reachable for a 6-aa segment.

### Independent cross-check

tcrdist3 ships its own hand-curated `db/adaptive_imgt_mapping.csv` (688 rows). It covers all 161
tokens and is **allele-identical on 155**. The 6 differences all favour the CDR-validated map:

| token | this table | tcrdist3 | why ours |
|---|---|---|---|
| `TCRAJ24-01` | `TRAJ24*02` | `TRAJ24*01` | the records end `…SWGKLQF` = `*02`; `*01` is `…SWGKFEF` |
| `TCRBV05-X` | `TRBV5-6` | `TRBV5-1` | observed `SGHDT`/`YYEEEE` = TRBV5-6 |
| `TCRBV06-X` | `TRBV6-5` | `TRBV6-1` | observed modal `MNHEY`/`SVGAGI` = TRBV6-5 |
| `TCRBV07-X` | `TRBV7-6` | `TRBV7-1` | observed `SGHVS`/`FNYEAQ` = TRBV7-6 |
| `TCRBV12-X` | `TRBV12-3` | `TRBV12-1` | observed `SGHNS`/`FNNNVP` = TRBV12-3 |
| `TCRBV12-03/12-04` | `TRBV12-3` | `TRBV12-1` | same |

tcrdist3 resolves a family call to the *first* member of the family; that is an arbitrary
tie-break, and for five of these six the first member is not the gene whose CDRs are actually there.

### Two defects in the source release, found by the same check

The CDR evidence also flags rows whose Adaptive gene call contradicts their own CDRs — these are
IMMREP25 annotation defects, recorded in the `evidence` column, not vdjtools bugs:

* `TCRBV06-05` (311 records): 190 rows carry `MNHEY`/`SVGAGI` (= `TRBV6-5`, consistent) but
  **120 rows carry `MNHNS`/`SASEGT`, which is `TRBV6-1`** — a different gene.
* `TCRBV06-X` (130 records): genuinely heterogeneous, as a family call should be — 70 rows TRBV6-5,
  40 rows TRBV6-2/6-3, 20 rows another TRBV6-5/6-6 allele.

## 4. Does the Java/Groovy original have the same defect?

**Yes, and slightly worse.** `legacy-1.x:src/main/groovy/…/misc/CommonUtil.groovy`:

```groovy
def res = it.replaceAll("TCR", "TR")
(1..9).each { int i -> res = res.replaceAll("0$i".toString(), "$i".toString()) }
```

Same global zero-strip, so the same subgroup/allele confusion, the same untouched `-X` and slash
ties, the same missing `/DVn`. Two extra hazards the Python port does not have: the patterns have no
digit boundary and are applied **sequentially** for `i = 1…9`, so an earlier rewrite can create a
match for a later one. The Python `re.sub(r"0([1-9])", r"\1", …)` is a single pass and is therefore
the strictly safer of the two — but both are wrong in the same way, so this is a **v1 bug inherited
by v2**, not a porting regression. Any published analysis that read immunoSEQ tables through legacy
vdjtools and then joined on gene name has the same silent losses.

## 5. Minimal implementation inside `convert.py`

Table lookup first, current behaviour as the fallback, so unknown tokens behave exactly as today.

```python
import csv
from functools import lru_cache
from importlib import resources


@lru_cache(maxsize=1)
def _adaptive_map() -> dict[str, str]:
    """Adaptive token → IMGT gene, from the CDR-validated shipped table."""
    txt = resources.files("vdjtools.resources").joinpath("adaptive_imgt_map.tsv").read_text()
    return {r["adaptive_token"]: r["imgt_gene"]
            for r in csv.DictReader(txt.splitlines(), delimiter="\t")
            if r["imgt_gene"] and r["status"] != "unresolved_dropped"}


def _adaptive_to_imgt(field: str | None) -> str | None:
    if field is None:
        return None
    tok = field.split(",")[0].replace('"', "").strip()      # first comma-tie, keep the slash tie
    if not tok or tok.lower() == "unresolved":
        return None
    hit = _adaptive_map().get(tok.split("*")[0])
    if hit:
        return hit
    gene = extract_vdj(field)                                # unchanged legacy path
    if gene is None or gene.lower() == "unresolved":
        return None
    return _ZERO_PAD.sub(r"\1", gene.replace("TCR", "TR"))
```

Note the one behavioural subtlety: the slash tie must be looked up **before** any splitting, because
`TCRBV03-01/03-02` is a key in the table. Keep `extract_vdj` untouched — it is shared by every other
reader. This body was executed verbatim against the shipped table (2026-08-08): all 10 self-check
assertions of §7 pass, plus `f("TCRBJ02-06*01") == "TRBJ2-6"` (allele suffix on the token),
`f("TCRBV99-99") == "TRBV99-99"` (off-table → legacy fallback, behaviour unchanged) and
`f(None) is None`.

Two follow-ups worth considering, neither required for the fix:

* Add a `strict: bool = False` keyword that returns `None` instead of falling back, so a caller can
  see how much of an immunoSEQ file is off-table rather than getting invented gene names.
* Longer term, resolve off-table tokens against arda's germline gene list at read time (the
  candidate-enumeration step of §3 needs nothing but the list of gene names) and keep the shipped
  table only for the CDR-arbitrated cases. This removes the table's coverage limit — it covers the
  161 tokens observed in IMMREP25 + pairSEQ, which is most but not all of Adaptive's vocabulary.
  Real immunoSEQ exports do **not** ship CDR1/CDR2 columns, so the CDR arbitration can only ever
  happen offline, in the builder — which is why a table is the right shape for this.

### Blast radius

* `convert._adaptive_to_imgt` — the function itself.
* `convert._adaptive_call(gene, family, ties)` — its only caller; the gene→family→ties fallback is
  unaffected, each level just resolves correctly now.
* `convert.read_immunoseq` — sets `V_CALL`, `D_CALL`, `J_CALL` from `_adaptive_call`. **All three**
  change, D included (`TCRBD01-01 → TRBD1`, not `TRBD1-1`).
* `io.batch` — `read` / `sniff_format` route `fmt="immunoseq"` here; `io.ingest_cohort(fmt=…)`
  inherits the change for every immunoSEQ cohort.
* Downstream of the gene call: anything joining on `v_call`/`j_call` — `model` (Pgen, germline
  resolution), `features` V/J usage, `overlap`, `preprocess.batch` V-J usage correction,
  `biomarker`. Corrected names will now *resolve* where they previously silently did not.
* **Not affected**: every other reader. `extract_vdj` is untouched; MiXcr/MiGEC/MiTCR/IMGT/Vidjil/
  RTCR/TRUST4/arda inputs already carry IMGT names.
* Existing tests stay green: `tests/python/test_convert.py` pins `immunoseq`/`immunoseqv2` to
  `v="TRBV29-1"`, `j="TRBJ2-6"` — both are `TCRBV29-01`/`TCRBJ02-06` in the fixture and both are
  unchanged by the table. Worth adding a fixture row with a `TCRAJ*-01` or `TCRBV09-01` call, since
  no current fixture exercises a failing token.
* **SOURCES.md** gets the table's provenance row (done).

## 6. Regenerating the table

```sh
# needs nothing but a stdlib python3; the two data inputs are read-only
python3 appendix/build_adaptive_imgt_map.py \
  --db /opt/homebrew/anaconda3/envs/cmp-tcrdist/lib/python3.8/site-packages/tcrdist/db/alphabeta_gammadelta_db.tsv \
  --immrep ~/vcs/manuscripts/2026-immrep25-audit/dump/immrep25/immrep2025_for_release.tsv \
  --tokens "$HOME/vcs/manuscripts/2026-immrep25-audit/results/pairseq_mock.tsv:va,ja,vb,jb" \
  --out python/vdjtools/resources/adaptive_imgt_map.tsv
```

Locate the tcrdist3 germline table with
`python -c "import tcrdist, os; print(os.path.dirname(tcrdist.__file__))"` in the `cmp-tcrdist`
conda env. Any IMGT-derived gene list with per-allele CDR1/CDR2 works in its place — the builder only
needs the columns `id`, `organism`, `chain`, `region`, `cdrs`.

## 7. Self-check (run after any change to `_adaptive_to_imgt`)

Assert-based, no test framework. **Three of the first five assertions fail on the current
implementation** — that is the point of the note.

```python
# save as appendix/check_adaptive_imgt.py, run from the repo root: python3 appendix/check_adaptive_imgt.py
import sys
sys.path.insert(0, "python")
from vdjtools.io.convert import _adaptive_to_imgt as f

# --- the five the fix must not regress ---------------------------------------------------
assert f("TCRBV29-01") == "TRBV29-1", f("TCRBV29-01")   # subgroup reading IS right here
assert f("TCRAV01-01") == "TRAV1-1", f("TCRAV01-01")    # ditto
assert f("TCRAJ39-01") == "TRAJ39", f("TCRAJ39-01")     # FAILS today -> "TRAJ39-1"
assert f("TCRAV22-01") == "TRAV22", f("TCRAV22-01")     # FAILS today -> "TRAV22-1"
assert f("unresolved") is None, f("unresolved")

# --- the other three failure modes -------------------------------------------------------
assert f("TCRBV09-01") == "TRBV9", f("TCRBV09-01")               # mode C, 795 calls
assert f("TCRBD01-01") == "TRBD1", f("TCRBD01-01")               # mode C, D genes
assert f("TCRBV03-01/03-02") == "TRBV3-1", f("TCRBV03-01/03-02")  # mode A, slash tie
assert f("TCRBV12-X") == "TRBV12-3", f("TCRBV12-X")              # mode B, family call
assert f("TCRAV38-02") == "TRAV38-2/DV8", f("TCRAV38-02")        # mode D, co-locus name

# --- the table must never emit a name outside the IMGT reference --------------------------
# (optional; needs the germline list)
# from vdjtools.model.reference import load_germline
# known = {a.split("*")[0] for a in load_germline("human")}
# for tok in ("TCRAJ39-01", "TCRBV09-01", "TCRBD01-01", "TCRAV38-02"):
#     assert f(tok) in known, (tok, f(tok))

print("adaptive->IMGT self-check OK")
```

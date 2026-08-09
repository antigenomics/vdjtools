#!/usr/bin/env python3
# 2026-08-08
# Build a CDR-VALIDATED Adaptive/immunoSEQ -> IMGT allele map.
#
# Why: vdjtools' `_adaptive_to_imgt` normalises Adaptive gene names with a *global* regex
# `0([1-9]) -> \1`, which silently invents genes that do not exist (TCRAJ39-01 -> "TRAJ39-1",
# TCRBV09-01 -> "TRBV9-1", TCRAV22-01 -> "TRAV22-1").  Whether the trailing "-01" of an
# Adaptive token is an IMGT *subgroup* or an *allele* is a per-gene-family fact, so no textual
# rule can decide it.  This script decides it from data: the IMGT gene list supplies the
# candidates, and the germline CDR1/CDR2 amino-acid strings that the IMMREP25 release ships
# alongside each record select the candidate that is actually right.
#
# Inputs (all plain text, stdlib only -- no third-party imports):
#   --db      tcrdist3 germline table `alphabeta_gammadelta_db.tsv` (columns id/organism/chain/
#             region/.../cdr_columns/cdrs; `cdrs` = "CDR1;CDR2;CDR2.5;CDR3-start", IMGT-gapped)
#   --immrep  IMMREP25 release TSV, supplies the CDR evidence
#             (tcra_v/tcra_j/tcra_cdr1/tcra_cdr2/tcra_cdr3 and the tcrb_* twins)
#   --tokens  optional extra TSVs contributing tokens but no CDRs (e.g. the pairSEQ mock);
#             given as PATH:va,ja,vb,jb  (column names for alpha-V, alpha-J, beta-V, beta-J)
#
# Output: one row per Adaptive token with the chosen IMGT allele, the evidence and a status.
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict

CDR_KEEP = str.maketrans("", "", ".-")


def ungap(s):
    return s.translate(CDR_KEEP)


# --- germline table -----------------------------------------------------------------------

def load_db(path, organism):
    """id -> dict(gene, allele, chain, region, cdr1, cdr2, cdr25, cdr3seg), ungapped."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["organism"] != organism:
                continue
            cdrs = (r.get("cdrs") or "").split(";")
            cdrs += [""] * (4 - len(cdrs))
            gene, _, allele = r["id"].partition("*")
            # V rows carry four gapped CDRs (CDR1;CDR2;CDR2.5;CDR3-start); J rows carry a single
            # field -- the J-side CDR3/junction segment through the conserved Phe/Trp118.
            out[r["id"]] = dict(gene=gene, allele=allele, chain=r["chain"], region=r["region"],
                               cdr1=ungap(cdrs[0]), cdr2=ungap(cdrs[1]),
                               cdr25=ungap(cdrs[2]),
                               cdr3seg=ungap(cdrs[0] if r["region"] == "J" else cdrs[3]))
    return out


_CORE = re.compile(r"^TR([ABDG])([VDJ])(\d+)(?:-(\d+))?$")


def gene_core(gene):
    """('TRAV14/DV4') -> ('TRAV14', locus='A', seg='V', fam=14, sub=None, orphon=False).

    The `/DVn` co-locus tag and the `/ORn-n` orphon tag are suffixes on an otherwise regular
    IMGT name; strip them to get a parseable core, but keep the orphon flag so orphons are
    never offered as candidates for a plain family call.
    """
    orphon = "/OR" in gene
    core = re.sub(r"/(DV\d+|OR\d+-\d+)$", "", gene)
    m = _CORE.match(core)
    if not m:
        return None
    loc, seg, fam, sub = m.groups()
    return dict(core=core, locus=loc, seg=seg, fam=int(fam),
                sub=None if sub is None else int(sub), orphon=orphon)


# --- Adaptive token parsing ---------------------------------------------------------------

_TOK = re.compile(r"^TCR([ABDG])([VDJ])(\d+)(?:-(\d+|X|or\d+_\d+))?$", re.I)


def parse_token(tok):
    """Adaptive token -> list of (locus, seg, fam, sub) calls; sub is int, 'X' or ('or',f,s).

    A slash-separated token ("TCRBV03-01/03-02") is an Adaptive *tie*: several equally-scored
    genes.  Every tie member becomes its own call and they are resolved jointly.
    """
    if not tok or not tok.upper().startswith("TCR"):
        return None
    parts = tok.split("/")
    head = _TOK.match(parts[0])
    if not head:
        return None
    loc, seg, fam, sub = head.groups()
    loc, seg = loc.upper(), seg.upper()

    def mksub(s):
        if s is None:
            return None
        if s.upper() == "X":
            return "X"
        m = re.match(r"or(\d+)_(\d+)$", s, re.I)
        if m:
            return ("or", int(m.group(1)), int(m.group(2)))
        return int(s)

    calls = [(loc, seg, int(fam), mksub(sub))]
    for p in parts[1:]:
        m = re.match(r"^(\d+)-(\d+|X)$", p, re.I)
        if m:
            calls.append((loc, seg, int(m.group(1)), mksub(m.group(2))))
        else:  # e.g. "TCRBV20-or09_02" already consumed; anything else is unparseable
            return None
    return calls


def candidate_genes(cores, call):
    """IMGT genes compatible with one Adaptive call.

    The trailing number of an Adaptive token is ambiguous by construction, so both readings
    are offered and the gene list / CDRs arbitrate:
      * subgroup reading  TCRAV01-01 -> TRAV1-1
      * allele reading    TCRAJ39-01 -> TRAJ39  (TRAJ has no subgroups at all)
    `-X` is Adaptive's "family only" call -> every gene of that family.
    `-orNN_MM` is an orphon call -> the /ORNN-MM gene.
    """
    loc, seg, fam, sub = call
    fam_genes = [g for g, c in cores.items() if c["locus"] == loc and c["seg"] == seg
                 and c["fam"] == fam]
    if sub == "X":
        return sorted(g for g in fam_genes if not cores[g]["orphon"])
    if isinstance(sub, tuple):
        want = "/OR%d-%d" % (sub[1], sub[2])
        return sorted(g for g in fam_genes if g.endswith(want))
    out = [g for g in fam_genes if not cores[g]["orphon"] and cores[g]["sub"] == sub]  # subgroup
    out += [g for g in fam_genes if not cores[g]["orphon"] and cores[g]["sub"] is None]  # allele
    return sorted(set(out))


# --- evidence -----------------------------------------------------------------------------

SLOTS = [("A", "V", "tcra_v", "tcra_cdr1", "tcra_cdr2", "tcra_cdr3"),
         ("A", "J", "tcra_j", None, None, "tcra_cdr3"),
         ("B", "V", "tcrb_v", "tcrb_cdr1", "tcrb_cdr2", "tcrb_cdr3"),
         ("B", "J", "tcrb_j", None, None, "tcrb_cdr3")]


def read_immrep(path):
    """(chain, seg, token) -> Counter of (cdr1, cdr2) and Counter of cdr3, plus record counts."""
    cdrs = defaultdict(Counter)
    cdr3 = defaultdict(Counter)
    nrec = Counter()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for chain, seg, tcol, c1, c2, c3 in SLOTS:
                tok = (r.get(tcol) or "").strip()
                if not tok:
                    continue
                k = (chain, seg, tok)
                nrec[k] += 1
                if c1:
                    cdrs[k][((r.get(c1) or "").strip(), (r.get(c2) or "").strip())] += 1
                cdr3[k][(r.get(c3) or "").strip()] += 1
    return cdrs, cdr3, nrec


def read_extra(spec):
    """PATH:va,ja,vb,jb -> ((chain,seg,token) -> n, (chain,seg,token) -> Counter(cdr3))."""
    path, _, cols = spec.rpartition(":")
    va, ja, vb, jb = cols.split(",")
    nrec, cdr3 = Counter(), defaultdict(Counter)
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        c3a = "cdr3a" if "cdr3a" in (rdr.fieldnames or []) else None
        c3b = "cdr3b" if "cdr3b" in (rdr.fieldnames or []) else None
        for r in rdr:
            for chain, seg, tcol, c3 in (("A", "V", va, None), ("A", "J", ja, c3a),
                                         ("B", "V", vb, None), ("B", "J", jb, c3b)):
                tok = (r.get(tcol) or "").strip()
                if tok:
                    nrec[(chain, seg, tok)] += 1
                    if c3:
                        cdr3[(chain, seg, tok)][(r.get(c3) or "").strip()] += 1
    return nrec, cdr3


# --- resolution ---------------------------------------------------------------------------

def resolve(token, chain, seg, db, cores, by_gene, obs_cdrs, obs_cdr3):
    """-> dict of output fields for one Adaptive token."""
    res = dict(chain=chain, segment=seg, adaptive_token=token, imgt_allele="", imgt_gene="",
               evidence="", n_candidates=0, cdr1_match="", cdr2_match="",
               allele_ambiguous="", status="no_candidate")
    if token.lower() in ("unresolved", "", "na", "nan"):
        res["evidence"] = "adaptive placeholder"
        res["status"] = "unresolved_dropped"
        return res
    calls = parse_token(token)
    if calls is None:
        res["evidence"] = "token does not parse"
        return res

    cands = sorted({g for c in calls for g in candidate_genes(cores, c)})
    alleles = sorted(a for g in cands for a in by_gene[g])
    res["n_candidates"] = len(alleles)
    if not alleles:
        res["evidence"] = "no gene of this family/subgroup in the IMGT reference"
        return res

    tie = len(calls) > 1 or (calls[0][3] == "X")

    # ---- V: match the germline CDR1/CDR2 the release ships for records carrying this token
    if seg == "V" and obs_cdrs:
        n_all = sum(obs_cdrs.values())
        top, n_top = obs_cdrs.most_common(1)[0]
        hit = [a for a in alleles if db[a]["cdr1"] == top[0] and db[a]["cdr2"] == top[1]]
        res["cdr1_match"] = str(any(db[a]["cdr1"] == top[0] for a in alleles))
        res["cdr2_match"] = str(any(db[a]["cdr2"] == top[1] for a in alleles))
        # Every *other* observed CDR pair is either another allele of the same gene (harmless
        # allelic variation) or another gene entirely (a heterogeneous token: an Adaptive `-X`
        # family call, an Adaptive tie, or an annotation error in the source release).
        other = []
        for pair, n in obs_cdrs.most_common()[1:]:
            gs = sorted({db[a]["gene"] for a in alleles
                         if db[a]["cdr1"] == pair[0] and db[a]["cdr2"] == pair[1]})
            if not gs:  # not any candidate -- name the gene it really is, reference-wide
                gs = sorted({d["gene"] for a, d in db.items() if d["region"] == "V"
                             and d["cdr1"] == pair[0] and d["cdr2"] == pair[1]})
                gs = ["!" + g for g in gs]  # "!" = outside the token's candidate set
            other.append((n, pair, gs))
        if hit:
            genes = sorted({db[a]["gene"] for a in hit})
            ev = ["cdr1=%s cdr2=%s matches %d/%d candidate allele(s) in %d/%d rows"
                  % (top[0], top[1], len(hit), len(alleles), n_top, n_all)]
            if len(genes) > 1:
                ev.append("CDR-identical genes " + ",".join(genes))
            het = [(n, gs) for n, _, gs in other if not gs or set(gs) - set(genes)]
            if het:
                ev.append("heterogeneous: " + "; ".join(
                    "%d rows carry CDRs of %s" % (n, ",".join(gs) or "no candidate")
                    for n, gs in het))
            allelic = sum(n for n, _, gs in other if gs and not set(gs) - set(genes))
            if allelic:
                ev.append("%d further rows carry another allele of the same gene" % allelic)
            res.update(imgt_allele=hit[0], imgt_gene=db[hit[0]]["gene"], evidence="; ".join(ev),
                       allele_ambiguous=str(len(hit) > 1 or len(genes) > 1 or bool(allelic)),
                       status="tie_resolved" if tie else "validated_by_cdr")
            return res
        # CDRs available but no candidate reproduces them -> do not invent a mapping
        res.update(evidence="cdr1=%s cdr2=%s (%d/%d rows) matches none of %s"
                            % (top[0], top[1], n_top, n_all, ",".join(alleles)),
                   status="unresolved_dropped")
        return res

    # ---- J: CDR1/CDR2 do not exist for J genes.  The Adaptive J family number is unique in
    # IMGT (TRAJ/TRBJ have at most one gene per family number), so the gene list is the primary
    # evidence; the CDR3 3' end is the cross-check.  Only the germline *3' end* of the J segment
    # survives recombination -- the 5' end is nibbled and overwritten by junctional insertions --
    # so score candidates by the length of the common suffix, not by a whole-segment match.
    if seg == "J" and obs_cdr3:
        n_all = sum(obs_cdr3.values())

        def suffix_score(seg_aa):
            if not seg_aa:
                return 0.0
            tot = 0
            for s, n in obs_cdr3.items():
                k = 0
                while k < min(len(s), len(seg_aa)) and s[-1 - k] == seg_aa[-1 - k]:
                    k += 1
                tot += k * n
            return tot / n_all

        scored = sorted(((suffix_score(db[a]["cdr3seg"]), a) for a in alleles),
                        key=lambda t: (-t[0], t[1]))
        best_score, best = scored[0]
        genes = sorted({db[a]["gene"] for a in alleles})
        best_genes = sorted({db[a]["gene"] for s, a in scored if s == best_score})
        ev = ("mean CDR3 3'-end agreement with germline J segment %s = %.2f aa over %d rows"
              % (db[best]["cdr3seg"], best_score, n_all))
        # >=5 exact 3' residues (the Phe/Trp118 anchor plus four) is the acceptance bar.
        if best_score >= 5.0 and len(best_genes) == 1:
            res.update(imgt_allele=best, imgt_gene=db[best]["gene"], evidence=ev,
                       allele_ambiguous=str(sum(1 for s, _ in scored if s == best_score) > 1),
                       status="tie_resolved" if tie else "validated_by_cdr")
            return res
        if len(genes) == 1:  # gene list already unambiguous; CDR3 check merely weak
            res.update(imgt_allele=alleles[0], imgt_gene=genes[0],
                       evidence="unique J gene of this family in the IMGT reference; " + ev,
                       allele_ambiguous=str(len(alleles) > 1), status="gene_list_unique")
            return res
        res.update(evidence="CDR3 3'-end cannot separate candidate J genes (%s); " % ",".join(genes)
                            + ev, status="unresolved_dropped")
        return res

    # ---- no sequence evidence at all: accept only if the gene list is unambiguous
    genes = sorted({db[a]["gene"] for a in alleles})
    if len(genes) == 1:
        res.update(imgt_allele=alleles[0], imgt_gene=genes[0],
                   evidence="unique gene in the IMGT reference; no CDR evidence in the inputs",
                   allele_ambiguous=str(len(alleles) > 1),
                   status="tie_resolved" if tie else "gene_list_unique")
        return res
    res.update(evidence="%d candidate genes (%s), no CDR evidence to choose"
                        % (len(genes), ",".join(genes)), status="unresolved_dropped")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", required=True, help="tcrdist3 alphabeta_gammadelta_db.tsv")
    ap.add_argument("--immrep", required=True, help="IMMREP25 release TSV (CDR evidence)")
    ap.add_argument("--tokens", action="append", default=[],
                    help="extra token source, PATH:va_col,ja_col,vb_col,jb_col")
    ap.add_argument("--also", default="",
                    help="comma-separated literal tokens to include with no record support "
                         "(gene-list resolution only) -- e.g. the D-gene calls, which "
                         "`read_immunoseq` maps but no paired cohort carries")
    ap.add_argument("--organism", default="human")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sep", default="\t")
    args = ap.parse_args()

    db = load_db(args.db, args.organism)
    by_gene = defaultdict(list)
    for a, d in db.items():
        by_gene[d["gene"]].append(a)
    for g in by_gene:
        by_gene[g].sort()
    cores = {g: c for g in by_gene if (c := gene_core(g))}

    obs_cdrs, obs_cdr3, nrec = read_immrep(args.immrep)
    for spec in args.tokens:
        n2, c32 = read_extra(spec)
        nrec.update(n2)
        for k, v in c32.items():
            obs_cdr3[k].update(v)
    for tok in filter(None, (t.strip() for t in args.also.split(","))):
        calls = parse_token(tok)
        if calls:
            nrec.setdefault(("A" if calls[0][0] == "A" else "B", calls[0][1], tok), 0)

    rows = []
    for (chain, seg, tok) in sorted(nrec):
        r = resolve(tok, chain, seg, db, cores, by_gene,
                    obs_cdrs.get((chain, seg, tok)), obs_cdr3.get((chain, seg, tok)))
        r["n_records"] = nrec[(chain, seg, tok)]
        rows.append(r)

    cols = ["chain", "segment", "adaptive_token", "n_records", "imgt_allele", "imgt_gene",
            "evidence", "n_candidates", "cdr1_match", "cdr2_match", "allele_ambiguous", "status"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter=args.sep)
        w.writeheader()
        for r in rows:
            w.writerow({c: r[c] for c in cols})

    st = Counter(r["status"] for r in rows)
    print("wrote %s (%d tokens)" % (args.out, len(rows)), file=sys.stderr)
    for k, v in st.most_common():
        print("  %-20s %4d tokens %7d records"
              % (k, v, sum(r["n_records"] for r in rows if r["status"] == k)), file=sys.stderr)


if __name__ == "__main__":
    main()

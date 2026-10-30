"""Explicit post-selection mass scaling; this is not a fitted selection model."""
import numpy as np

# ALICE's published average on the OLGA TRB scale; not established for TRA.
ALICE_TRB_SELECTION_Q = 9.41


def selection_mass(pgen, selection_q=1.0):
    """Scale generation mass by an explicit positive selection factor.

    Returns Q * Pgen without clipping or renormalization. The result is the
    mass used in a generative expectation, not a calibrated binding probability.
    Empirical control frequencies already reflect selection; do not rescale them.
    """
    mass = np.asarray(pgen, dtype=float)
    q = float(selection_q)
    if not np.isfinite(q) or q <= 0 or not np.all(np.isfinite(mass) & (mass >= 0)):
        raise ValueError("finite nonnegative generation mass and positive selection_q required")
    with np.errstate(over="ignore"):
        scaled = mass * q
    if not np.all(np.isfinite(scaled)):
        raise ValueError("selection-scaled mass overflow")
    return scaled.item() if scaled.ndim == 0 else scaled


def fit_selection_factor(counts, pgen, *, total_events):
    """Fit Q to observed sequence multiplicities by zero-truncated Poisson MLE.

    Counts are independent receptor-event counts, never read abundances. Only
    observed sequences enter, so the likelihood conditions on count > 0.
    The unconditioned mean is total_events * Q * Pgen for each sequence.
    Zero-model-mass observations must be audited by the caller, not hidden here.
    """
    from scipy.optimize import brentq

    y = np.asarray(counts, dtype=float)
    p = np.asarray(pgen, dtype=float)
    n = float(total_events)
    if y.ndim != 1 or not len(y) or p.shape != y.shape:
        raise ValueError("nonempty aligned one-dimensional counts and pgen required")
    if not np.all(np.isfinite(y) & (y >= 1) & (y == np.floor(y))):
        raise ValueError("positive integer receptor-event counts required")
    if not np.isfinite(n) or n >= 2**53 or n < y.sum() or n != np.floor(n):
        raise ValueError("integer total_events below 2**53 must cover the observed counts")
    if not np.all(np.isfinite(p) & (p > 0) & (p <= 1)):
        raise ValueError("positive finite Pgen required; audit unsupported sequences separately")
    if y.sum() == len(y):
        raise ValueError("all counts are singletons; no positive selection-factor estimate")
    exposure = n * p
    def score(logq):
        lam = np.exp(logq) * exposure
        mean = np.divide(lam, -np.expm1(-lam), out=np.ones_like(lam), where=lam > 0)
        return float(mean.sum() - y.sum())
    lo, hi = -1., 1.
    while score(lo) > 0:
        lo *= 2
        if lo < -700: raise ValueError("selection factor below numerical range")
    while score(hi) < 0:
        hi *= 2
        if hi > 700: raise ValueError("selection factor above numerical range")
    return float(np.exp(brentq(score, lo, hi)))


def control_event_counts(path, *, sequence_column="cdr3aa", weight_column="incidence",
                         max_events=None, seed=42):
    """Stream a nucleotide-key control pool into amino-acid event counts.

    The caller establishes independent donors. NT incidence counts one event per
    donor/NT key; summing NT variants gives amino-acid occurrences, never reads.
    """
    import polars as pl
    from ..io import iter_table_batches

    if max_events is not None and (isinstance(max_events, bool) or max_events <= 0 or int(max_events) != max_events):
        raise ValueError("max_events must be a positive integer")
    parts = []
    rows = invalid_rows = events = invalid_events = 0
    for batch in iter_table_batches(path):
        columns = [weight_column] + [c for c in ("occurrences", "incidence", "convergence")
                                     if c in batch.columns and c != weight_column]
        for column in columns:
            text = batch[column].cast(pl.String)
            if not text.str.contains(r"^[0-9]+$").fill_null(False).all():
                raise ValueError(f"{column} requires nonnegative integer event counts")
        batch = batch.with_columns(pl.col(columns).cast(pl.Int64))
        if "convergence" in batch.columns and (batch["convergence"] != 1).any():
            raise ValueError("control input must be nucleotide-keyed, with convergence=1")
        if {"occurrences", "incidence"} <= set(batch.columns) and (batch["occurrences"] != batch["incidence"]).any():
            raise ValueError("NT occurrence/incidence differ; resolve event provenance before fitting")
        valid = pl.col(sequence_column).str.contains(r"^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$").fill_null(False)
        rejected = batch.filter(~valid)
        rows += batch.height
        events += int(batch[weight_column].sum())
        invalid_rows += rejected.height
        invalid_events += int(rejected[weight_column].sum())
        parts.append(batch.filter(valid).select(pl.col(sequence_column).alias("junction_aa"),
                     pl.col(weight_column).alias("events")).group_by("junction_aa").agg(pl.col("events").sum()))
    if not parts: raise ValueError("empty control pool")
    counts = pl.concat(parts).group_by("junction_aa").agg(pl.col("events").sum()).filter(pl.col("events") > 0).sort("junction_aa")
    if events >= 2**53 or int(counts["events"].sum()) != events-invalid_events:
        raise ValueError("event total exceeds supported precision or aggregation overflow")
    population_events = events - invalid_events
    if max_events is not None and max_events < population_events:
        # Exact event sampling without expanding incidence counts into individual rows.
        sampled = np.random.default_rng(seed).multivariate_hypergeometric(
            counts["events"].to_numpy(), int(max_events), method="marginals")
        counts = counts.with_columns(pl.Series("events", sampled)).filter(pl.col("events") > 0)
    return counts, dict(source_rows=rows, source_events=events, invalid_rows=invalid_rows,
                        invalid_events=invalid_events, retained_events=int(counts["events"].sum()),
                        population_events=population_events, max_events=max_events,
                        sampling_seed=seed if max_events is not None else None,
                        sampling="without replacement across receptor events" if max_events is not None else "full pool",
                        observed_junctions=counts.height, weight_column=weight_column)


def selection_factor(locus, *, organism="human", source="isalgo"):
    """Return the calibrated effective multiplier for an unchanged bundled model.

    Unsupported model/species/locus contexts require an explicit user factor.
    Model hashes are checked before returning Q; factors never transfer silently
    between parameter sets. Empirical control frequencies do not need this factor.
    """
    import hashlib
    import json
    from importlib.resources import files
    from .bundled import _source_family, _bundled_key

    family=_source_family(source)
    species=organism.lower()
    locus=locus.upper()
    root=files("vdjtools.model")
    records=json.loads((root/"selection_factors.json").read_text())
    record=next((r for r in records if (r["source"],r["species"],r["locus"]) ==
                 (family,species,locus)),None)
    if record is None:
        raise ValueError(f"no calibrated selection factor for {family}/{species}/{locus}; supply explicit Q")
    model=root/"_bundled"/family/_bundled_key(family,locus,species)
    actual={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in model.iterdir() if p.is_file()}
    if actual!=record["model_sha256"]:
        raise ValueError("bundled model changed since selection calibration; supply explicit Q or recalibrate")
    return selection_mass(1.,record["selection_q"])

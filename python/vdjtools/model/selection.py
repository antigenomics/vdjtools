"""Model-bound global selection calibration and optional coarse selection shapes."""
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
    q = np.asarray(selection_q, dtype=float)
    if q.ndim and q.shape != mass.shape:
        raise ValueError('selection_q must be scalar or aligned with generation mass')
    if not np.all(np.isfinite(q) & (q > 0)) or not np.all(np.isfinite(mass) & (mass >= 0)):
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
        if lo < -700:
            raise ValueError("selection factor below numerical range")
    while score(hi) < 0:
        hi *= 2
        if hi > 700:
            raise ValueError("selection factor above numerical range")
    return float(np.exp(brentq(score, lo, hi)))


def control_event_counts(path, *, sequence_column="cdr3aa", weight_column="incidence",
                         max_events=None, seed=42, preserve_vj=False):
    """Stream a nucleotide-key control pool into amino-acid event counts.

    The caller establishes independent donors. NT incidence counts one event per
    donor/NT key; summing NT variants gives amino-acid occurrences, never reads.
    """
    import polars as pl
    from ..io import iter_table_batches

    if max_events is not None and (isinstance(max_events, bool) or max_events <= 0 or int(max_events) != max_events):
        raise ValueError("max_events must be a positive integer")
    def validated(batch):
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
        if int(batch[weight_column].cast(pl.Int128).sum()) >= 2**53:
            raise ValueError('event total exceeds supported precision')
        return batch

    if preserve_vj:
        # Two raw passes bound retained gene strings by the event sample size.
        from polars import selectors as cs
        names = list(dict.fromkeys([sequence_column, weight_column, 'occurrences',
            'incidence', 'convergence', 'v_call', 'j_call', 'v', 'j']))
        def batches():
            pending = None
            for batch in iter_table_batches(path, columns=cs.by_name(names, require_all=False)):
                pending = batch if pending is None else pl.concat([pending, batch], rechunk=False)
                while pending.height >= 65536:
                    yield validated(pending.head(65536))
                    pending = pending.slice(65536)
            if pending is not None and pending.height:
                yield validated(pending)
        rows = invalid_rows = events = invalid_events = 0
        valid = pl.col(sequence_column).str.contains(r'^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$').fill_null(False)
        for batch in batches():
            rejected = batch.filter(~valid)
            rows += batch.height
            events += int(batch[weight_column].sum())
            invalid_rows += rejected.height
            invalid_events += int(rejected[weight_column].sum())
            if events >= 2**53:
                raise ValueError('event total exceeds supported precision')
        if not rows:
            raise ValueError('empty control pool')
        population = remaining = events - invalid_events
        if max_events is not None and max_events < population and population >= 10**9:
            raise ValueError('event subsampling requires population below 1000000000')
        take = min(population, int(max_events)) if max_events is not None else population
        rng = np.random.default_rng(seed)
        parts = []
        for batch in batches():
            batch = batch.filter(valid)
            weights = batch[weight_column].to_numpy()
            size = int(weights.sum())
            selected = size if take == remaining else (take if size == remaining else int(rng.hypergeometric(size, remaining-size, take)))
            if selected:
                sampled = weights if selected == size else rng.multivariate_hypergeometric(weights, selected, method='marginals')
                genes = [pl.col(c if c in batch.columns else c[0]).cast(pl.String).alias(c)
                         if c in batch.columns or c[0] in batch.columns else pl.lit(None, dtype=pl.String).alias(c)
                         for c in ('v_call', 'j_call')]
                parts.append(batch.with_columns(pl.Series('_events', sampled)).filter(pl.col('_events') > 0)
                    .select(pl.col(sequence_column).alias('junction_aa'), *genes, pl.col('_events').alias('events'))
                    .group_by('junction_aa', 'v_call', 'j_call').agg(pl.col('events').sum()))
            remaining -= size
            take -= selected
        if remaining or take:
            raise RuntimeError('control event population changed while sampling')
        counts = (pl.concat(parts).group_by('junction_aa', 'v_call', 'j_call').agg(pl.col('events').sum())
                  .sort('junction_aa', 'v_call', 'j_call')) if parts else pl.DataFrame(
                      schema={'junction_aa': pl.String, 'v_call': pl.String, 'j_call': pl.String, 'events': pl.Int64})
        return counts, dict(source_rows=rows, source_events=events, invalid_rows=invalid_rows,
            invalid_events=invalid_events, population_events=population, retained_events=int(counts['events'].sum()),
            max_events=max_events, sampling_seed=seed if max_events is not None else None,
            sampling='streaming-hypergeometric-v1' if max_events is not None else 'full pool',
            observed_junctions=counts['junction_aa'].n_unique(), observed_strata=counts.height, weight_column=weight_column)
    parts = []
    rows = invalid_rows = events = invalid_events = 0
    for batch in iter_table_batches(path):
        batch = validated(batch)
        valid = pl.col(sequence_column).str.contains(r"^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$").fill_null(False)
        rejected = batch.filter(~valid)
        rows += batch.height
        events += int(batch[weight_column].sum())
        invalid_rows += rejected.height
        invalid_events += int(rejected[weight_column].sum())
        parts.append(batch.filter(valid).select(pl.col(sequence_column).alias("junction_aa"),
                     pl.col(weight_column).alias("events")).group_by("junction_aa").agg(pl.col("events").sum()))
    if not parts:
        raise ValueError("empty control pool")
    counts = pl.concat(parts).group_by("junction_aa").agg(pl.col("events").sum()).filter(pl.col("events") > 0).sort("junction_aa")
    if events >= 2**53 or int(counts["events"].sum()) != events-invalid_events:
        raise ValueError("event total exceeds supported precision or aggregation overflow")
    population_events = events - invalid_events
    if max_events is not None and max_events < population_events:
        if population_events >= 10**9:
            raise ValueError('event subsampling requires population below 1000000000')
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


def selection_reference(model, *, threads=0):
    """Exact unnormalized eligible-junction mass by V gene, J gene and AA length.

    Eligibility is C-first, F/W-last and no stop codons, matching control-event
    ingestion. Length includes both anchors. Allele histories are summed into
    gene cells; cells with zero model mass are omitted. The sum of ``p_pre`` is
    the model's eligible probability, not necessarily one. Callers normalize
    by that sum when fitting a selection shape and keep its global scale separate.

    Enumerate the complete finite length support from packed germline/insertion
    bounds and score all cells in one native degenerate-AA batch. This computes
    no amino-acid-position selection parameters and uses no sampled background.
    """
    import polars as pl
    from ..cores import kernel_threads, require_memory
    from . import native

    packed, vi, ji = native.pack(model)
    vs = [name for name, i in vi.items() if i in packed.func_v]
    js = [name for name, i in ji.items() if i in packed.func_j]
    max_nt = max(map(len, packed.cut_v), default=0) + max(map(len, packed.cut_j), default=0)
    if packed.vdj:
        max_nt += max(map(len, packed.cut_d), default=0)
        max_nt += max(0, len(packed.ins_vd) - 1) + max(0, len(packed.ins_dj) - 1)
        if packed.dd and packed.p_nd2 > 0:
            max_nt += max(map(len, packed.cut_d), default=0) + max(0, len(packed.ins_dd) - 1)
    else:
        max_nt += max(0, len(packed.ins_vj) - 1)
    lengths = range(2, max_nt // 3 + 1)
    budget = kernel_threads(threads)
    n = len(vs) * len(js) * len(lengths)
    # Python residue sets, C++ strings/masks, flat output columns and live DP.
    require_memory(len(vs) * len(js) * sum(lengths) * 192 + n * 256 +
                   budget * (len(js) + len(packed.func_d) + 4) * (max_nt + 1) * 25 * 8,
                   'exact V/J/length selection reference')
    motifs = {length: ['C'] + ['ACDEFGHIKLMNPQRSTVWY'] * (length - 2) + ['FW'] for length in lengths}
    result = (pl.DataFrame({'v_call': vs}, schema={'v_call': pl.String})
              .join(pl.DataFrame({'j_call': js}, schema={'j_call': pl.String}), how='cross')
              .join(pl.DataFrame({'junction_aa_length': list(lengths)},
                                 schema={'junction_aa_length': pl.Int64}), how='cross'))
    probabilities = native.pgen_aa_degenerate_batch(model,
        [motifs[length] for length in result['junction_aa_length']],
        v=result['v_call'].to_list(), j=result['j_call'].to_list(), threads=budget)
    if not np.all(np.isfinite(probabilities) & (np.asarray(probabilities) >= 0)):
        raise RuntimeError('native Pgen returned invalid selection-reference mass')
    return (result.with_columns(pl.Series('p_pre', probabilities, dtype=pl.Float64),
                              pl.col('v_call', 'j_call').str.split('*').list.first())
            .group_by('v_call', 'j_call', 'junction_aa_length').agg(pl.col('p_pre').sum())
            .filter(pl.col('p_pre') > 0).sort('v_call', 'j_call', 'junction_aa_length'))


def fit_selection_profile(model, observations, *, kind='vj-length', smoothing=1.,
                          shrinkage=1., threads=0):
    """Fit a normalized coarse selection shape with a shared smooth length curve.

    ``observations`` contains integer ``junction_aa_length`` and ``events``;
    optional V/J calls may be comma-separated sets. Missing/unknown calls are
    marginalized, and incompatible preliminary calls are reconsidered at their
    observed length. Only lengths with no eligible model mass are unsupported.
    The penalized likelihood sums probabilities over each compatible gene set.
    ``smoothing`` penalizes second differences of log length weights;
    positive ``shrinkage`` penalizes squared log weights and fixes their gauges.
    Both penalties use event-count likelihood units. The returned shape has
    mean one under the exact eligible preselection reference; it has no global Q.
    """
    import polars as pl
    from scipy.optimize import minimize
    from scipy.special import logsumexp
    from ..cores import kernel_threads, numerical_kernels, require_memory
    from . import native
    from .pgen import _resolve_gene_call

    if kind not in {'length', 'vj-length'}:
        raise ValueError('profile kind must be length or vj-length')
    if not np.isfinite(smoothing) or smoothing < 0 or not np.isfinite(shrinkage) or shrinkage <= 0:
        raise ValueError('finite nonnegative smoothing and positive shrinkage required')
    if not {'junction_aa_length', 'events'} <= set(observations.columns) or not observations.height:
        raise ValueError('nonempty length/event observations required')
    for name in ('junction_aa_length', 'events'):
        if not observations[name].dtype.is_integer() or observations[name].null_count() or observations[name].min() < 1:
            raise ValueError(f'positive integer {name} required')
    if float(observations['events'].cast(pl.Float64).sum()) >= 2**53:
        raise ValueError('event total must be below 2**53')
    input_events = int(observations['events'].sum())
    budget = kernel_threads(threads)
    reference = selection_reference(model, threads=budget)
    if not reference.height:
        raise ValueError('model has no eligible selection-reference mass')
    eligible_mass = float(reference['p_pre'].sum())
    obs = observations.select('junction_aa_length', 'events', *[
        pl.col(name).cast(pl.String) if name in observations.columns else pl.lit(None, dtype=pl.String).alias(name)
        for name in ('v_call', 'j_call')])
    _, vi, ji = native.pack(model)
    for name, index in (('v_call', vi), ('j_call', ji)):
        known = set(reference[name])
        mapping = {}
        for raw in obs[name].drop_nulls().unique():
            members = {_resolve_gene_call(index, part.strip()).split('*')[0]
                       for part in raw.split(',') if part.strip()}
            mapping[raw] = ','.join(sorted(members)) if members and members <= known else None
        obs = obs.with_columns(pl.col(name).replace_strict(mapping, default=None, return_dtype=pl.String))
    unresolved = int(obs.filter(pl.col('v_call').is_null() | pl.col('j_call').is_null())['events'].sum())
    if kind == 'length':
        fit_reference = reference.group_by('junction_aa_length').agg(pl.col('p_pre').sum()).with_columns(
            pl.lit('*').alias('v_call'), pl.lit('*').alias('j_call'))
        obs = obs.with_columns(pl.lit(None, dtype=pl.String).alias('v_call'),
                               pl.lit(None, dtype=pl.String).alias('j_call'))
    else:
        fit_reference = reference
    obs = obs.group_by('junction_aa_length', 'v_call', 'j_call').agg(pl.col('events').sum()).with_row_index('_obs')
    lo = int(fit_reference['junction_aa_length'].min())
    hi = int(fit_reference['junction_aa_length'].max())
    pairs = fit_reference.select('v_call', 'j_call').unique().sort('v_call', 'j_call').with_row_index('_gene')
    cells = fit_reference.with_row_index('_cell').join(pairs, on=['v_call', 'j_call'],
                                                      validate='m:1', maintain_order='left')
    expansion = int(obs.select(((pl.col('v_call').str.count_matches(',').fill_null(0) + 1) *
                               (pl.col('j_call').str.count_matches(',').fill_null(0) + 1)).sum()).item())
    require_memory(expansion * 160 + cells.height * 160, 'selection gene-set expansion')
    expanded = (obs.with_columns(pl.col('v_call', 'j_call').str.split(','))
                .explode('v_call', empty_as_null=True).explode('j_call', empty_as_null=True))
    edges = []
    edge_total = 0
    for has_v, has_j in ((True, True), (True, False), (False, True), (False, False)):
        part = expanded.filter((pl.col('v_call').is_not_null() == has_v) &
                               (pl.col('j_call').is_not_null() == has_j))
        keys = ['junction_aa_length'] + (['v_call'] if has_v else []) + (['j_call'] if has_j else [])
        sizes = part.join(cells.group_by(keys).len('_cells'), on=keys, how='inner')
        n_edges = int(sizes['_cells'].sum())
        edge_total += n_edges
        require_memory(edge_total * 160 + cells.height * 160 + obs.height * 128,
                       'selection observation likelihood')
        edges.append(part.select('_obs', *keys).join(cells, on=keys, how='inner').select('_obs', '_cell'))
    edges = pl.concat(edges).unique()
    missing = obs.join(edges.select('_obs').unique(), on='_obs', how='anti')
    reassigned = missing.join(cells.select('junction_aa_length').unique(), on='junction_aa_length', how='semi')
    fallback_size = int(reassigned.join(cells.group_by('junction_aa_length').len('_cells'),
        on='junction_aa_length', how='inner')['_cells'].sum())
    require_memory((edge_total + fallback_size) * 160 + cells.height * 160 + obs.height * 128,
                   'selection reassigned observation likelihood')
    edges = pl.concat([edges, reassigned.select('_obs', 'junction_aa_length')
                       .join(cells, on='junction_aa_length').select('_obs', '_cell')]).sort('_obs', '_cell')
    supported = edges['_obs'].unique().sort().to_numpy()
    if not len(supported):
        raise ValueError(f'no supported observations among {input_events} events')
    remap = np.full(obs.height, -1, dtype=np.intp)
    remap[supported] = np.arange(len(supported))
    oi = remap[edges['_obs'].to_numpy()]
    ci = edges['_cell'].to_numpy().astype(np.intp)
    y = obs['events'].to_numpy()[supported].astype(float)
    n = float(y.sum())
    gi = cells['_gene'].to_numpy().astype(np.intp)
    li = cells['junction_aa_length'].to_numpy().astype(np.intp) - lo
    p0 = cells['p_pre'].to_numpy() / eligible_mass
    logp0 = np.log(p0)
    ng = pairs.height if kind == 'vj-length' else 0
    nl = hi - lo + 1

    def objective(theta):
        a = theta[:ng] if ng else np.zeros(1)
        b = theta[ng:]
        logcell = logp0 + a[gi] + b[li]
        logz = logsumexp(logcell)
        maxima = np.full(len(y), -np.inf)
        np.maximum.at(maxima, oi, logcell[ci])
        edge = np.exp(logcell[ci] - maxima[oi])
        sums = np.bincount(oi, weights=edge, minlength=len(y))
        logden = maxima + np.log(sums)
        posterior = np.bincount(ci, weights=y[oi] * edge / sums[oi], minlength=len(p0))
        residual = n * np.exp(logcell - logz) - posterior
        gradient = np.r_[np.bincount(gi, weights=residual, minlength=ng) if ng else [],
                         np.bincount(li, weights=residual, minlength=nl)]
        second = np.diff(b, n=2)
        penalty = .5 * (shrinkage * np.dot(theta, theta) + smoothing * np.dot(second, second))
        gradient += shrinkage * theta
        gradient[ng:ng+nl-2] += smoothing * second
        gradient[ng+1:ng+nl-1] -= 2 * smoothing * second
        gradient[ng+2:] += smoothing * second
        return (n * logz - np.dot(y, logden) + penalty) / n, gradient / n

    with numerical_kernels(budget):
        result = minimize(objective, np.zeros(ng + nl), jac=True, method='L-BFGS-B',
                          options={'maxiter': 2000, 'ftol': 1e-12, 'gtol': 1e-9})
    if not result.success or not np.all(np.isfinite(result.x)):
        raise RuntimeError(f'selection-profile optimization failed: {result.message}')
    a = result.x[:ng] if ng else np.zeros(1)
    b = result.x[ng:]
    profile = dict(version=1, kind=kind, parameter_sha256=_selection_model_hash(model), length_min=lo, log_length=b.tolist(),
        log_vj=[dict(v_call=v, j_call=j, log_q=float(a[i])) for i, v, j in pairs.iter_rows()] if ng else [],
        log_normalizer=float(logsumexp(logp0 + a[gi] + b[li])), eligible_mass=eligible_mass,
        smoothing=float(smoothing), shrinkage=float(shrinkage),
        fit=dict(iterations=int(result.nit), penalized_nll_per_event=float(result.fun),
                 max_gradient_per_event=float(np.max(np.abs(result.jac)))),
        audit=dict(input_events=input_events, supported_events=int(n), unsupported_events=input_events-int(n),
                   input_strata=obs.height, supported_strata=len(supported),
                   unresolved_gene_events=unresolved, reassigned_gene_events=int(reassigned['events'].sum())))
    factors = selection_profile_factors(profile, reference['junction_aa_length'].to_numpy(),
        reference['v_call'].to_list(), reference['j_call'].to_list())
    profile['selected_vj'] = (reference.with_columns(
        pl.Series('p', reference['p_pre'].to_numpy() * factors / eligible_mass))
        .group_by('v_call', 'j_call').agg(pl.col('p').sum()).sort('v_call', 'j_call').to_dicts())
    return profile


def selection_profile_factors(profile, lengths, v=None, j=None):
    """Evaluate a fitted shape on canonical V/J genes and anchored AA lengths."""
    import polars as pl

    values = np.asarray(lengths, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values) & (values == np.floor(values))):
        raise ValueError('one-dimensional integer AA lengths required')
    curve = np.asarray(profile['log_length'], dtype=float)
    index = values.astype(np.intp) - int(profile['length_min'])
    if np.any((index < 0) | (index >= len(curve))):
        raise ValueError('length outside fitted selection-profile support')
    logweight = curve[index] - float(profile['log_normalizer'])
    if profile['kind'] == 'vj-length':
        if v is None or j is None or len(v) != len(values) or len(j) != len(values):
            raise ValueError('aligned V/J gene calls required for vj-length factors')
        queries = pl.DataFrame({'v_call': v, 'j_call': j},
                               schema={'v_call': pl.String, 'j_call': pl.String}).with_columns(
            pl.col('v_call', 'j_call').str.split('*').list.first())
        matched = queries.join(pl.DataFrame(profile['log_vj']), on=['v_call', 'j_call'],
                               how='left', validate='m:1', maintain_order='left')
        if matched['log_q'].null_count():
            raise ValueError('V/J pair outside fitted selection-profile support')
        logweight = logweight + matched['log_q'].to_numpy()
    elif profile['kind'] != 'length':
        raise ValueError('profile kind must be length or vj-length')
    factors = np.exp(logweight)
    if not np.all(np.isfinite(factors) & (factors > 0)):
        raise ValueError('selection factor outside floating-point range')
    return factors


def _selection_model_hash(model):
    """Portable parameter identity, using the existing corpus raw-table hash convention."""
    import hashlib
    import json
    from dataclasses import asdict
    digest = hashlib.sha256(json.dumps(dict(locus=model.locus, organism=model.organism,
        chain_type=model.chain_type, palindrome_max=model.manifest.palindrome_max,
        events={name: asdict(event) for name, event in model.manifest.events.items()}), sort_keys=True).encode())
    for collection in (model.tables, model.genomic):
        for name, table in sorted(collection.items()):
            columns = sorted(table.columns)
            digest.update((name + '\n').encode())
            digest.update(table.select(columns).sort(columns).write_csv(
                float_scientific=True, float_precision=17).encode())
    return digest.hexdigest()


def selection_pgen_aa_batch(model, sequences, profile, *, threads=0):
    """Exact marginal generation mass times a fitted selection shape, without global Q.

    Sum gene-dependent weights inside the latent V/J sum. The model is checked
    against the profile's parameter identity and remains unchanged. This is an
    explicitly selected mass, not the original model's Pgen.
    """
    import polars as pl
    from ..cores import kernel_threads, require_memory
    from . import native

    if profile['parameter_sha256'] != _selection_model_hash(model):
        raise ValueError('model parameters differ from selection-profile calibration')
    if not len(sequences):
        return []
    seq = pl.Series('junction_aa', sequences, dtype=pl.String).str.to_uppercase()
    if not seq.str.contains(r'^C[ACDEFGHIKLMNPQRSTVWY]*[FW]$').fill_null(False).all():
        raise ValueError('selection profiles require anchored stop-free AA junctions')
    lengths = seq.str.len_chars().to_numpy()
    curve = np.asarray(profile['log_length'])
    li = lengths.astype(np.intp) - profile['length_min']
    if np.any((li < 0) | (li >= len(curve))):
        raise ValueError('length outside fitted selection-profile support')
    length_weights = np.exp(curve[li] - profile['log_normalizer'])
    if not np.all(np.isfinite(length_weights) & (length_weights > 0)):
        raise ValueError('selection length factor outside floating-point range')
    budget = kernel_threads(threads)
    if profile['kind'] == 'length':
        values = native.pgen_aa_batch(model, seq.to_list(), threads=budget)
    elif profile['kind'] == 'vj-length':
        from .._core import selection_pgen_aa_batch as _selected
        packed, vi, ji = native.pack(model)
        parameter_bytes = sum(table.estimated_size() for table in [*model.tables.values(), *model.genomic.values()])
        require_memory(len(ji) * parameter_bytes * 4 + len(seq) * (int(lengths.max()) * 16 + 256),
                       'selection-weighted native Pgen')
        v = pl.DataFrame({'v_call': list(vi), '_v': list(vi.values())})
        j = pl.DataFrame({'j_call': list(ji), '_j': list(ji.values())})
        weights = (v.with_columns(pl.col('v_call').str.split('*').list.first())
                   .join(j.with_columns(pl.col('j_call').str.split('*').list.first()), how='cross')
                   .join(pl.DataFrame(profile['log_vj']), on=['v_call', 'j_call'], how='left', validate='m:1')
                   .sort('_v', '_j').select(pl.col('log_q').exp().fill_null(0).alias('weight'))['weight'].to_list())
        values = _selected(packed, seq.to_list(), weights, budget)
    else:
        raise ValueError('profile kind must be length or vj-length')
    result = np.asarray(values) * length_weights
    if not np.all(np.isfinite(result) & (result >= 0)):
        raise ValueError('invalid selection-weighted generation mass')
    return result.tolist()


def load_selection_profile(report, model):
    """Load an explicit calibration report and verify its complete model identity."""
    import json
    from pathlib import Path
    identity = _selection_model_hash(model)
    if isinstance(report, str) and report in {'length', 'vj-length'}:
        from importlib.resources import files
        records = json.loads((files('vdjtools.model') / 'selection_profiles.json').read_text())
        record = next((r for r in records if r['selection_profile']['kind'] == report and
                       r['selection_profile']['parameter_sha256'] == identity), None)
        if record is None:
            raise ValueError(f'no bundled {report} selection profile for this model; supply a calibration report')
        report = record
    elif not isinstance(report, dict):
        report = json.loads(Path(report).read_text())
    profile = report['selection_profile']
    if profile.get('version') != 1 or profile.get('kind') not in {'length', 'vj-length'}:
        raise ValueError('unsupported selection-profile schema')
    if profile['parameter_sha256'] != identity:
        raise ValueError('model parameters differ from selection-profile calibration')
    q = selection_mass(1., report['selection_q'])
    return profile, q

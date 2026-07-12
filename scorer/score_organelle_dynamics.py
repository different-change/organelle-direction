"""
score_organelle_dynamics — directional organelle-dynamics scoring from expression.

A species-agnostic, single-sample module scorer that separates organelle
BIOGENESIS (build-up / proliferation) from organelle DEGRADATION (selective
autophagy / clearance), and returns a signed net-direction score per sample.

Motivation
----------
A generic "peroxisome activity" signature cannot tell proliferation apart from
pexophagy — both move many of the same genes. This function scores two directional
modules independently and reports their difference, so the *sign* of the net score
distinguishes a shift toward biogenesis (positive) from a shift toward clearance
(negative). It also reports a DEGRADATION-SELECTIVE score built only from the
organelle-selective degradation genes (e.g. the pexophagy receptor ATG36), because
the core autophagy machinery overlaps bulk autophagy and dilutes specificity.

Design choices
--------------
- Single-sample, background-normalized scoring. Default method='rank' ranks all
  genes *within each sample* and reports the mean percentile-rank of module genes
  (0..1). This needs no cross-sample information, so a score is portable to a lone
  sample and comparable across platforms. method='zscore' (z across samples, mean z
  of module genes) is provided for the validation comparison where all samples of a
  study are scored together.
- Empirical null. For each module, `n_perm` random gene sets of the same size are
  drawn from the background (all genes present in the matrix) to give an empirical
  z-score and two-sided p per sample — so a module score is interpretable versus
  chance rather than on an arbitrary scale.
- Robust to missing genes; coverage (module genes found / requested) is reported.

mRNA caveat
-----------
Transcript abundance is a proxy. A positive net_direction is "consistent with a
shift toward biogenesis"; a negative one "consistent with a shift toward clearance".
These are NOT organelle turnover rates.

Author: Phase 1, directional peroxisome-dynamics project.
"""

from __future__ import annotations
import numpy as np
import pandas as pd

__all__ = ["score_organelle_dynamics"]


# ----------------------------------------------------------------------------- #
# helpers
# ----------------------------------------------------------------------------- #
def _looks_like_log(expr: pd.DataFrame) -> bool:
    """Heuristic: is the matrix already on a log scale?"""
    finite = expr.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return True
    # Log-scale data rarely exceeds ~40; raw counts/intensities routinely reach 1e3-1e6.
    return np.nanmax(finite) < 40 and np.nanmin(finite) > -40


def _prepare_matrix(expr: pd.DataFrame, id_map, assume_log) -> pd.DataFrame:
    """Map IDs to the module namespace, coerce numeric, log-normalize defensively."""
    m = expr.copy()
    if id_map is not None:
        m.index = m.index.map(lambda x: id_map.get(x, x))
    # collapse duplicate IDs (e.g. multiple probes) by median
    m = m.apply(pd.to_numeric, errors="coerce")
    if m.index.duplicated().any():
        m = m.groupby(level=0).median()
    # drop all-NaN rows
    m = m.dropna(how="all")
    is_log = assume_log if assume_log is not None else _looks_like_log(m)
    if not is_log:
        # log2(x+1) after clipping negatives (defensive); counts/intensities only
        m = np.log2(m.clip(lower=0) + 1.0)
    return m, is_log


def _percentile_ranks(col: pd.Series) -> pd.Series:
    """Within-sample percentile rank in 0..1 (average ties)."""
    return col.rank(method="average", pct=True)


def _module_stat(ranks_or_z: pd.Series, module_ids) -> float:
    present = [g for g in module_ids if g in ranks_or_z.index]
    if not present:
        return np.nan
    return float(ranks_or_z.loc[present].mean())


def _empirical_null(sample_vec: pd.Series, k: int, n_perm: int, rng) -> np.ndarray:
    """Draw n_perm random gene-set means of size k from the background."""
    bg = sample_vec.dropna().to_numpy()
    if k == 0 or bg.size < k:
        return np.array([])
    idx = rng.integers(0, bg.size, size=(n_perm, k))
    return bg[idx].mean(axis=1)


# ----------------------------------------------------------------------------- #
# main
# ----------------------------------------------------------------------------- #
def score_organelle_dynamics(
    expr: pd.DataFrame,
    biogenesis_ids,
    degradation_ids,
    *,
    degradation_selective_ids=None,
    id_map=None,
    method: str = "rank",
    n_perm: int = 1000,
    random_state: int = 0,
    assume_log=None,
) -> pd.DataFrame:
    """
    Score directional organelle dynamics per sample.

    Parameters
    ----------
    expr : DataFrame (genes x samples)
        Index = gene IDs in the *module* namespace, or supply `id_map` to translate
        the matrix index into that namespace. Values may be raw or log; the function
        log2(x+1)-normalizes if the matrix does not already look log-scaled (override
        with `assume_log`).
    biogenesis_ids, degradation_ids : iterable of gene IDs
        The two directional modules (module namespace).
    degradation_selective_ids : iterable, optional
        Organelle-SELECTIVE degradation genes (e.g. pexophagy receptor ATG36) — the
        subset of `degradation_ids` not shared with bulk autophagy. Used for the
        `degradation_selective_score`. If None, this column is NaN.
    id_map : dict, optional
        Maps matrix index -> module namespace (applied before scoring).
    method : {'rank','zscore'}
        'rank' (default): within-sample mean percentile-rank of module genes (0..1);
        single-sample, portable. 'zscore': gene z across samples, mean z of module
        genes; needs >=2 samples.
    n_perm : int
        Number of random background gene-sets for the empirical null (per module,
        per sample). 0 disables the null (empirical z / p columns become NaN).
    random_state : int
        Seed for the permutation null.
    assume_log : bool, optional
        Force the log-scale assumption instead of auto-detecting.

    Returns
    -------
    DataFrame indexed by sample, with columns:
        biogenesis_score, degradation_score, degradation_selective_score,
        net_direction (= biogenesis - degradation),
        net_direction_selective (= biogenesis - degradation_selective),
        <module>_emp_z, <module>_emp_p  (empirical z / two-sided p vs null),
        biogenesis_coverage, degradation_coverage, degradation_selective_coverage
        (found / requested), n_background (genes used as the null pool).

    Notes
    -----
    Pure function: no module lists are hard-coded here; pass them in. mRNA is a
    proxy — interpret sign as "consistent with a shift toward biogenesis/clearance",
    never as an organelle turnover rate.
    """
    if method not in ("rank", "zscore"):
        raise ValueError("method must be 'rank' or 'zscore'")

    bio = list(dict.fromkeys(biogenesis_ids))
    deg = list(dict.fromkeys(degradation_ids))
    deg_sel = list(dict.fromkeys(degradation_selective_ids)) if degradation_selective_ids is not None else None

    m, is_log = _prepare_matrix(expr, id_map, assume_log)
    if m.shape[1] == 0:
        raise ValueError("expression matrix has no samples")
    if method == "zscore" and m.shape[1] < 2:
        raise ValueError("method='zscore' needs >=2 samples; use method='rank' for single samples")

    rng = np.random.default_rng(random_state)
    samples = list(m.columns)

    # coverage (relative to the matrix's available genes)
    present_idx = set(m.index)
    cov = {
        "biogenesis": (len([g for g in bio if g in present_idx]), len(bio)),
        "degradation": (len([g for g in deg if g in present_idx]), len(deg)),
        "degradation_selective": (
            (len([g for g in deg_sel if g in present_idx]), len(deg_sel)) if deg_sel is not None else (0, 0)
        ),
    }

    # per-sample transformed vector used for scoring + null
    if method == "rank":
        transformed = {s: _percentile_ranks(m[s]) for s in samples}
    else:  # zscore: z per gene across samples
        mu = m.mean(axis=1)
        sd = m.std(axis=1, ddof=0).replace(0, np.nan)
        z = m.sub(mu, axis=0).div(sd, axis=0)
        transformed = {s: z[s] for s in samples}

    rows = []
    for s in samples:
        vec = transformed[s]
        rec = {}
        for mod_name, ids in [("biogenesis", bio), ("degradation", deg),
                              ("degradation_selective", deg_sel if deg_sel is not None else [])]:
            present = [g for g in ids if g in vec.index and np.isfinite(vec.get(g, np.nan))]
            stat = float(vec.loc[present].mean()) if present else np.nan
            rec[f"{mod_name}_score"] = stat
            # empirical null
            if n_perm and present:
                null = _empirical_null(vec, len(present), n_perm, rng)
                if null.size:
                    nmu, nsd = null.mean(), null.std(ddof=0)
                    ez = (stat - nmu) / nsd if nsd > 0 else np.nan
                    # two-sided empirical p
                    p = (np.sum(np.abs(null - nmu) >= abs(stat - nmu)) + 1) / (null.size + 1)
                    rec[f"{mod_name}_emp_z"] = ez
                    rec[f"{mod_name}_emp_p"] = p
                else:
                    rec[f"{mod_name}_emp_z"] = np.nan
                    rec[f"{mod_name}_emp_p"] = np.nan
            else:
                rec[f"{mod_name}_emp_z"] = np.nan
                rec[f"{mod_name}_emp_p"] = np.nan
        rec["net_direction"] = rec["biogenesis_score"] - rec["degradation_score"]
        rec["net_direction_selective"] = rec["biogenesis_score"] - rec["degradation_selective_score"]
        rec["n_background"] = int(vec.dropna().shape[0])
        rows.append(rec)

    out = pd.DataFrame(rows, index=samples)
    for k, (found, req) in cov.items():
        out[f"{k}_coverage"] = f"{found}/{req}"
        out[f"{k}_coverage_frac"] = (found / req) if req else np.nan
    out.attrs["method"] = method
    out.attrs["log_scale_used"] = is_log
    # tidy column order
    lead = ["biogenesis_score", "degradation_score", "degradation_selective_score",
            "net_direction", "net_direction_selective"]
    rest = [c for c in out.columns if c not in lead]
    return out[lead + rest]


# ----------------------------------------------------------------------------- #
# smoke test
# ----------------------------------------------------------------------------- #
if __name__ == "__main__":
    rng = np.random.default_rng(0)
    genes = [f"G{i:04d}" for i in range(2000)]
    bio_ids = genes[:20]
    deg_ids = genes[20:40]
    deg_sel = genes[20:23]
    # 3 "biogenesis-high" samples, 3 "degradation-high" samples
    base = pd.DataFrame(rng.normal(8, 1, size=(2000, 6)),
                        index=genes, columns=[f"bioUp{i}" for i in range(3)] + [f"degUp{i}" for i in range(3)])
    base.loc[bio_ids, [c for c in base.columns if c.startswith("bioUp")]] += 4
    base.loc[deg_ids, [c for c in base.columns if c.startswith("degUp")]] += 4
    res = score_organelle_dynamics(base, bio_ids, deg_ids,
                                   degradation_selective_ids=deg_sel, n_perm=500, assume_log=True)
    print(res[["biogenesis_score", "degradation_score", "net_direction",
               "biogenesis_emp_z", "degradation_emp_z"]].round(3))
    assert res.loc[[c for c in base.columns if c.startswith("bioUp")], "net_direction"].mean() > 0
    assert res.loc[[c for c in base.columns if c.startswith("degUp")], "net_direction"].mean() < 0
    print("\nsmoke test PASSED: biogenesis-up samples net>0, degradation-up samples net<0")

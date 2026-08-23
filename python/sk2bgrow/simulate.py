"""Monte-Carlo harness reproducing the design report's section 5.

The report's central quantitative claim is that a deterministic 16-enzyme anchor
union beats an equal-density random sketch, and beats a single enzyme badly, at
every coverage level — and that in the wet-lab route the union wins even at a
*fixed read budget*, because averaging within a window quenches site-efficiency
noise at rate sqrt(n).

This module re-runs those simulations so the claim is reproducible rather than
quoted. The estimator used for the comparison is deliberately the *Pilea-style*
pipeline (25 kb windows -> positive counts -> ZTP rate -> sort -> trimmed
log-linear fit), so that the anchor set is the only thing that varies. The
sk2bGrow V-shape estimator is available through ``estimator="v_shape"`` to show
separately what the new estimator adds on top of the new sketch.

Anchor coordinates come from the Rust layer (``sk2bgrow digest``) via
:func:`anchors_from_digest`, or from :func:`synthetic_anchors` when no genome is
at hand. The enzyme logic is never reimplemented here — one source of truth.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .fit import circular_distance, fit_v_shape
from .ztp import solve_ztp_lambda

__all__ = [
    "SimConfig",
    "replication_log2_copy",
    "simulate_counts",
    "estimate_sorted",
    "estimate",
    "synthetic_anchors",
    "anchors_from_digest",
    "run_grid",
    "route_a",
    "route_b",
]

#: E. coli K-12 MG1655, the genome the report's simulations use.
ECOLI_LEN = 4_641_652
ECOLI_ORI = 3_923_883

#: Anchor counts measured on E. coli K-12 in report section 4.1.
ECOLI_UNION_ANCHORS = 28_381
ECOLI_BCGI_ANCHORS = 2_932
#: Pilea's default sketch size on a genome this size (FracMinHash s=250).
PILEA_SKETCH_ANCHORS = 18_600


@dataclass
class SimConfig:
    """One simulated condition."""

    log2_ptr: float = 1.0
    genome_len: int = ECOLI_LEN
    ori: int = ECOLI_ORI
    #: Expected reads per anchor before the copy-number gradient is applied.
    per_site_depth: float = 1.0
    #: Lognormal sd of per-site efficiency (digestion / ligation / amplification).
    sigma_eff: float = 0.15
    window_bp: int = 25_000
    trim: float = 0.05
    n_reps: int = 150
    seed: int = 0


def replication_log2_copy(positions: np.ndarray, ori: float, genome_len: float, log2_ptr: float) -> np.ndarray:
    """Log2 relative copy number along the chromosome.

    The standard model: bidirectional replication makes copy number fall linearly
    in log space from the origin to the terminus, by ``log2_ptr`` in total.
    """
    d = circular_distance(positions, ori, genome_len)
    return -log2_ptr * (d / (genome_len / 2.0))


def simulate_counts(
    positions: np.ndarray,
    cfg: SimConfig,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw Poisson counts at each anchor.

    ``mean = per_site_depth * copy_number * efficiency``, with efficiency
    lognormal with unit mean (the ``-sigma^2/2`` term keeps the mean at 1, so
    changing ``sigma_eff`` changes the noise without also changing the depth).
    """
    log2_copy = replication_log2_copy(positions, cfg.ori, cfg.genome_len, cfg.log2_ptr)
    eff = (
        rng.lognormal(-0.5 * cfg.sigma_eff**2, cfg.sigma_eff, size=positions.size)
        if cfg.sigma_eff > 0
        else np.ones(positions.size)
    )
    lam = cfg.per_site_depth * np.exp2(log2_copy) * eff
    return rng.poisson(lam)


def estimate_sorted(positions: np.ndarray, counts: np.ndarray, cfg: SimConfig) -> float:
    """The report's reference pipeline, faithful to iRep/Pilea.

    Per 25 kb window: take the positive counts, invert the zero-truncated Poisson
    mean to recover the window rate, then sort the log2 rates, fit a line against
    rank with both ends trimmed, and multiply the slope by the window count.
    """
    bins = positions // cfg.window_bp
    rates = []
    for b in np.unique(bins):
        c = counts[bins == b]
        pos = c[c >= 1]
        if pos.size < 3:
            continue
        lam = solve_ztp_lambda(float(pos.mean()))
        if lam > 0:
            rates.append(lam)
    if len(rates) < 10:
        return np.nan
    y = np.sort(np.log2(np.array(rates)))
    n = y.size
    lo = int(np.floor(n * cfg.trim))
    hi = n - lo
    ys = y[lo:hi]
    if ys.size < 4:
        return np.nan
    slope = float(np.polyfit(np.arange(ys.size, dtype=float), ys, 1)[0])
    return slope * n


def estimate(positions: np.ndarray, counts: np.ndarray, cfg: SimConfig, estimator: str = "sorted") -> float:
    """Dispatch to the parity estimator or the sk2bGrow coordinate fit."""
    if estimator == "sorted":
        return estimate_sorted(positions, counts, cfg)
    if estimator != "v_shape":
        raise ValueError(f"unknown estimator {estimator!r}")

    bins = positions // cfg.window_bp
    mids, rates, ses = [], [], []
    for b in np.unique(bins):
        sel = bins == b
        c = counts[sel]
        pos = c[c >= 1]
        if pos.size < 3:
            continue
        lam = solve_ztp_lambda(float(pos.mean()))
        if lam <= 0:
            continue
        mids.append(float(positions[sel].mean()))
        rates.append(np.log2(lam))
        # Poisson-scale error propagated to log2.
        ses.append(1.0 / (np.sqrt(max(lam, 1e-9) * pos.size) * np.log(2.0)))
    if len(rates) < 10:
        return np.nan
    fit = fit_v_shape(np.array(mids), np.array(rates), cfg.genome_len, se=np.array(ses))
    return fit.log2_ptr


def synthetic_anchors(n: int, genome_len: int, kind: str = "uniform", seed: int = 0) -> np.ndarray:
    """Anchor coordinates for a sketch of ``n`` loci.

    ``kind="uniform"`` places them by Bernoulli sampling — the FracMinHash
    behaviour, whose spacings are geometric with no upper bound. ``kind="even"``
    spaces them regularly, which is the idealised deterministic limit. Real
    enzyme anchors sit between the two, which is why the report measures them on
    a real genome rather than assuming either.
    """
    rng = np.random.default_rng(seed)
    if kind == "uniform":
        return np.sort(rng.integers(0, genome_len, size=n))
    if kind == "even":
        return np.sort((np.arange(n) * (genome_len / n)).astype(np.int64))
    raise ValueError(f"unknown anchor kind {kind!r}")


def anchors_from_digest(path: str | Path, enzyme: str | None = None) -> np.ndarray:
    """Read anchor coordinates from a TGT text dump written by ``sk2bgrow index --write-tgt``.

    Passing ``enzyme`` restricts to one enzyme's anchors, which is how the
    single-enzyme arm of the report's comparison is built from real coordinates.
    """
    positions = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            name = fields[0].split(":", 1)[0]
            if enzyme is not None and name != enzyme:
                continue
            positions.append(int(fields[2]))
    return np.sort(np.array(positions, dtype=np.int64))


def run_grid(
    anchor_sets: dict[str, np.ndarray],
    depths: dict[str, float] | list[float],
    log2_ptrs: list[float],
    base: SimConfig | None = None,
    estimator: str = "sorted",
    sigma_eff: float | dict[str, float] | None = None,
) -> pd.DataFrame:
    """Run every (anchor set x depth x true PTR) cell and summarise bias and RMSE.

    ``depths`` may be a per-anchor-set dict, which is how the fixed-read-budget
    comparison is expressed: one enzyme at 50x per site and sixteen enzymes at
    the same total budget spread thinner.
    """
    base = base or SimConfig()
    rows = []
    for set_name, positions in anchor_sets.items():
        depth_values = depths[set_name] if isinstance(depths, dict) else depths
        if np.isscalar(depth_values):
            depth_values = [float(depth_values)]
        sig = sigma_eff.get(set_name, base.sigma_eff) if isinstance(sigma_eff, dict) else (
            base.sigma_eff if sigma_eff is None else float(sigma_eff)
        )
        for depth in np.atleast_1d(depth_values):
            for truth in log2_ptrs:
                cfg = SimConfig(
                    log2_ptr=truth,
                    genome_len=base.genome_len,
                    ori=base.ori,
                    per_site_depth=float(depth),
                    sigma_eff=sig,
                    window_bp=base.window_bp,
                    trim=base.trim,
                    n_reps=base.n_reps,
                    seed=base.seed,
                )
                ests = np.array(
                    [
                        estimate(positions, simulate_counts(positions, cfg, np.random.default_rng(cfg.seed + 1000 * r)), cfg, estimator)
                        for r in range(cfg.n_reps)
                    ]
                )
                ok = np.isfinite(ests)
                err = ests[ok] - truth
                rows.append(
                    {
                        "anchor_set": set_name,
                        "n_anchors": int(positions.size),
                        "per_site_depth": float(depth),
                        "sigma_eff": float(sig),
                        "true_log2_ptr": truth,
                        "estimator": estimator,
                        "n_ok": int(ok.sum()),
                        "n_reps": cfg.n_reps,
                        "bias": float(np.mean(err)) if err.size else np.nan,
                        "rmse": float(np.sqrt(np.mean(err**2))) if err.size else np.nan,
                        "sd": float(np.std(err)) if err.size else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def route_a(
    coverages: list[float] | None = None,
    log2_ptrs: list[float] | None = None,
    n_reps: int = 150,
    estimator: str = "sorted",
    seed: int = 0,
    anchor_sets: dict[str, np.ndarray] | None = None,
) -> pd.DataFrame:
    """Route A (WMS in silico): single enzyme vs random sketch vs 16-enzyme union.

    Per-site depth is ``0.8 x`` per-base coverage — the fraction of a 150 bp read
    that can contain a ~30 bp tag, which is the same 0.8 factor Pilea's k=31
    sketch carries. Comparing at matched *coverage* rather than matched
    per-anchor depth is what makes this an honest head-to-head.
    """
    coverages = coverages or [0.5, 1.0, 2.0, 5.0, 10.0]
    log2_ptrs = log2_ptrs or [1.0, 2.0]
    sets = anchor_sets or {
        "BcgI_single": synthetic_anchors(ECOLI_BCGI_ANCHORS, ECOLI_LEN, "uniform", seed + 1),
        "random_sketch": synthetic_anchors(PILEA_SKETCH_ANCHORS, ECOLI_LEN, "uniform", seed + 2),
        "union_16": synthetic_anchors(ECOLI_UNION_ANCHORS, ECOLI_LEN, "even", seed + 3),
    }
    base = SimConfig(sigma_eff=0.15, n_reps=n_reps, seed=seed)
    depths = {name: [0.8 * c for c in coverages] for name in sets}
    out = run_grid(sets, depths, log2_ptrs, base=base, estimator=estimator)
    out["coverage"] = out["per_site_depth"] / 0.8
    return out


def route_b(
    sigma_effs: list[float] | None = None,
    log2_ptrs: list[float] | None = None,
    single_depth: float = 50.0,
    n_reps: int = 150,
    estimator: str = "sorted",
    seed: int = 0,
) -> pd.DataFrame:
    """Route B (wet-lab 2bRAD): one enzyme deep vs sixteen enzymes at the same budget.

    The single-enzyme arm gets ``single_depth`` reads per site. The union arm
    gets the *same total reads* spread over ~9.7x as many anchors, so its
    per-site depth falls to about 5.2x. The 16x-budget arm shows what a full
    sixteen-library experiment would buy.
    """
    sigma_effs = sigma_effs or [0.3, 0.6]
    log2_ptrs = log2_ptrs or [1.0]
    single = synthetic_anchors(ECOLI_BCGI_ANCHORS, ECOLI_LEN, "uniform", seed + 1)
    union = synthetic_anchors(ECOLI_UNION_ANCHORS, ECOLI_LEN, "even", seed + 3)
    shared_budget = single_depth * ECOLI_BCGI_ANCHORS / ECOLI_UNION_ANCHORS

    frames = []
    for sig in sigma_effs:
        sets = {"BcgI_single": single, "union_16_same_budget": union, "union_16_full_budget": union}
        depths = {
            "BcgI_single": single_depth,
            "union_16_same_budget": shared_budget,
            "union_16_full_budget": single_depth,
        }
        base = SimConfig(sigma_eff=sig, n_reps=n_reps, seed=seed)
        frames.append(run_grid(sets, depths, log2_ptrs, base=base, estimator=estimator, sigma_eff=sig))
    return pd.concat(frames, ignore_index=True)

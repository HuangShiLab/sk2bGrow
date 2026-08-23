"""Multi-sample PTR dynamics.

The design report's §3.1 observation is that a *deterministic* sketch changes the
shape of a longitudinal analysis. With a random sketch, every sample is a
separately-drawn set of loci, so a time series is "estimate a noisy curve per
sample, then compare curves". With 2bRAD anchors the loci are physically fixed by
the enzymes, so the same genomic positions are re-observed in every sample and the
data become a repeated-measures table on fixed loci — the object that
mixed-effects models were built for, and the same paradigm digital karyotyping
established for copy number two decades ago.

Two views are provided:

* :func:`delta_ptr` — the summary view: genome x sample log2(PTR), differenced
  against a baseline with propagated errors.
* :func:`anchor_matrix` — the raw view: anchor x sample counts, which is what a
  DEMIC-style Poisson-PCA check or a repeated-measures model consumes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

from . import io as sk_io

__all__ = ["read_outputs", "ptr_matrix", "delta_ptr", "anchor_matrix", "trend_test"]


def read_outputs(paths: Iterable[str | Path]) -> pd.DataFrame:
    """Concatenate several ``output.tsv`` files."""
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t", na_values=["NA"])
        if "log2(PTR)" not in df.columns:
            raise sk_io.Sk2bIoError(f"{p} is not an sk2bgrow output.tsv (no 'log2(PTR)' column)")
        df["source"] = str(p)
        frames.append(df)
    if not frames:
        raise sk_io.Sk2bIoError("no output files given")
    return pd.concat(frames, ignore_index=True)


def ptr_matrix(outputs: pd.DataFrame, value: str = "log2(PTR)", qc_only: bool = True) -> pd.DataFrame:
    """Pivot to a genome x sample matrix."""
    df = outputs
    if qc_only and "pass_qc" in df.columns:
        df = df[df["pass_qc"].astype(str).str.lower().isin(["true", "1"])]
    return df.pivot_table(index="genome", columns="sample", values=value, aggfunc="mean")


def delta_ptr(
    outputs: pd.DataFrame,
    baseline: str | None = None,
    metadata: pd.DataFrame | None = None,
    qc_only: bool = True,
) -> pd.DataFrame:
    """Difference each sample's log2(PTR) against a baseline, per genome.

    ``baseline`` names a sample, or a group in ``metadata`` (columns
    ``sample``, ``group``, optional ``timepoint``). Without one, the per-genome
    mean across samples is the reference, which answers "which samples deviate"
    rather than "how much did it change from time zero".

    Errors are propagated as ``sqrt(se_i^2 + se_baseline^2)``: the two estimates
    are independent samples, so the variances add. The z-test that follows is a
    normal approximation, appropriate because each log2(PTR) is itself an
    inverse-variance-weighted mean over enzymes and windows.
    """
    df = outputs
    if qc_only and "pass_qc" in df.columns:
        df = df[df["pass_qc"].astype(str).str.lower().isin(["true", "1"])]
    df = df.copy()
    if "se" not in df.columns:
        df["se"] = np.nan

    if metadata is not None:
        df = df.merge(metadata, on="sample", how="left")

    rows = []
    for genome, grp in df.groupby("genome", sort=True):
        if baseline is not None and "group" in grp.columns and (grp["group"] == baseline).any():
            base = grp[grp["group"] == baseline]
        elif baseline is not None and (grp["sample"] == baseline).any():
            base = grp[grp["sample"] == baseline]
        elif baseline is not None:
            continue  # this genome was not measured in the baseline
        else:
            base = grp

        b_val = float(np.nanmean(base["log2(PTR)"]))
        b_se_vals = base["se"].to_numpy(dtype=float)
        b_se = float(np.sqrt(np.nansum(b_se_vals**2)) / max(np.isfinite(b_se_vals).sum(), 1))

        for _, r in grp.iterrows():
            if baseline is not None and r["sample"] in set(base["sample"]):
                continue
            d = float(r["log2(PTR)"]) - b_val
            se = float(np.sqrt(np.nan_to_num(r["se"], nan=0.0) ** 2 + b_se**2))
            z = d / se if se > 0 else np.nan
            rows.append(
                {
                    "genome": genome,
                    "sample": r["sample"],
                    "group": r.get("group", "-"),
                    "timepoint": r.get("timepoint", np.nan),
                    "log2_ptr": float(r["log2(PTR)"]),
                    "baseline_log2_ptr": b_val,
                    "delta_log2_ptr": d,
                    "se": se,
                    "z": z,
                    "p": float(2.0 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                    "fold_change": float(2.0**d),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["genome", "sample"]).reset_index(drop=True)
        out["q"] = _bh(out["p"].to_numpy(dtype=float))
    return out


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR. NaNs pass through untouched."""
    q = np.full_like(p, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    n = ranked.size
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def anchor_matrix(count_tables: Iterable[str | Path], genome: str | None = None, usable_only: bool = True) -> pd.DataFrame:
    """Build the anchor x sample count matrix.

    Rows are anchors keyed by ``genome:contig:position:enzyme`` — stable across
    samples because the enzymes, not a random hash, chose them. That stability is
    the whole reason this matrix is well defined.
    """
    df = sk_io.concat_counts(count_tables)
    if genome is not None:
        df = df[df["genome"] == genome]
    if usable_only:
        df = df[df["usable"]]
    if df.empty:
        return pd.DataFrame()
    df = df.copy()
    df["anchor"] = (
        df["genome"].astype(str)
        + ":"
        + df["contig_id"].astype(str)
        + ":"
        + df["position"].astype(str)
        + ":"
        + df["enzyme"].astype(str)
    )
    return df.pivot_table(index="anchor", columns="sample", values="count", aggfunc="sum", fill_value=0)


def trend_test(deltas: pd.DataFrame, time_column: str = "timepoint") -> pd.DataFrame:
    """Per-genome linear trend of log2(PTR) against time.

    A weighted least-squares slope with inverse-variance weights, plus its
    standard error and a t-test. This is the simplest form of the repeated-
    measures analysis the fixed-anchor design enables; a full mixed-effects model
    with an anchor random effect belongs on top of
    :func:`anchor_matrix`.
    """
    rows = []
    for genome, grp in deltas.groupby("genome", sort=True):
        g = grp[np.isfinite(grp[time_column].astype(float))] if time_column in grp else grp.iloc[0:0]
        if len(g) < 3:
            rows.append({"genome": genome, "slope": np.nan, "se": np.nan, "p": np.nan, "n": len(g),
                         "note": "fewer than 3 timepoints"})
            continue
        t = g[time_column].to_numpy(dtype=float)
        y = g["log2_ptr"].to_numpy(dtype=float)
        s = g["se"].to_numpy(dtype=float)
        w = np.where(np.isfinite(s) & (s > 0), 1.0 / s**2, np.nan)
        w = np.where(np.isfinite(w), w, np.nanmedian(w) if np.isfinite(np.nanmedian(w)) else 1.0)
        design = np.column_stack([np.ones_like(t), t])
        sw = np.sqrt(w)
        beta, *_ = np.linalg.lstsq(design * sw[:, None], y * sw, rcond=None)
        resid = y - design @ beta
        dof = max(len(t) - 2, 1)
        red = float(np.sum(w * resid**2) / dof)
        cov = np.linalg.pinv((design * sw[:, None]).T @ (design * sw[:, None])) * red
        se = float(np.sqrt(max(cov[1, 1], 0.0)))
        tstat = beta[1] / se if se > 0 else np.nan
        rows.append(
            {
                "genome": genome,
                "slope": float(beta[1]),
                "se": se,
                "t": float(tstat) if np.isfinite(tstat) else np.nan,
                "p": float(2.0 * stats.t.sf(abs(tstat), dof)) if np.isfinite(tstat) else np.nan,
                "n": int(len(t)),
                "note": "",
            }
        )
    return pd.DataFrame(rows)

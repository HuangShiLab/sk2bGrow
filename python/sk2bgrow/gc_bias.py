"""Per-enzyme GC bias correction.

Pilea applies one global loess correction to window coverage and drops samples
that stay overdispersed afterwards — report defect D6. Two things change here:

* **the curve is fitted inside each enzyme.** Every Type IIB recognition site has
  its own GC composition, so the 16 enzymes sample 16 different, individually
  narrow GC neighbourhoods. One global curve averages them into a shape that fits
  none of them.
* **the correction is applied at anchor resolution.** GC is stored per anchor
  (+/-250 bp), so the offset is computed per anchor and only then averaged into
  the window. Correcting on a window's mean GC would throw away exactly the
  within-window variation the correction exists to remove.

The correction is an *offset in log2 space*, never a rescaling of the counts:
the window models in :mod:`sk2bgrow.ztp` need integer counts, and a multiplicative
fudge would silently break the zero-truncation.

A note on what this does **not** need to fix: a constant per-enzyme efficiency
factor cannot bias PTR at all, because each enzyme is fitted separately and a
constant offset is absorbed by that fit's intercept. Only GC *slope* within an
enzyme matters. The per-enzyme factors are still reported, for QC.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess

__all__ = ["GcCurve", "fit_curves", "add_anchor_offsets", "apply_to_windows", "tukey_mask", "enzyme_efficiency"]

_PSEUDO = 0.5


@dataclass
class GcCurve:
    """A monotone-interpolated loess curve of log2 efficiency against GC."""

    enzyme: str
    gc: np.ndarray
    log2_offset: np.ndarray
    n_anchors: int
    #: Peak-to-trough amplitude of the fitted curve, in log2 units. A large
    #: amplitude means this enzyme's anchors are strongly GC-biased in this
    #: sample and its PTR estimate deserves less weight.
    amplitude: float
    #: Fraction of the raw curve that survived shrinkage, in [0, 1]. 0 means the
    #: curve explained no more variance than fitting noise would, and no
    #: correction is applied.
    shrinkage: float = 1.0
    #: Variance of log2 count explained by the curve.
    r2: float = np.nan

    def __call__(self, gc: np.ndarray) -> np.ndarray:
        """Evaluate the curve, clamping outside the fitted GC range.

        Extrapolating a loess fit is how a handful of extreme-GC anchors acquire
        wild corrections, so the ends are held flat instead.
        """
        gc = np.asarray(gc, dtype=float)
        out = np.interp(gc, self.gc, self.log2_offset, left=self.log2_offset[0], right=self.log2_offset[-1])
        return np.where(np.isfinite(gc), out, 0.0)


def fit_curves(counts: pd.DataFrame, frac: float = 0.4, min_anchors: int = 100, usable_only: bool = True) -> dict[str, GcCurve]:
    """Fit one loess curve per enzyme.

    Anchors with zero counts are kept: a GC neighbourhood that is systematically
    unobserved is exactly the bias being measured, and dropping zeros would hide
    it.

    Enzymes with fewer than ``min_anchors`` usable anchors get no curve — a loess
    through a handful of points fits noise, and no correction is safer than a
    wrong one.

    **Shrinkage.** A loess will happily trace Poisson scatter, and a sparse enzyme
    has few anchors to average over: fitting one on a genome with no GC bias at
    all produced spurious "corrections" of ~0.5 log2 in testing. Since each
    enzyme gets its own curve, and each curve then distorts its enzyme
    differently, that noise turns into *between-enzyme* disagreement and the
    cross-enzyme consistency test fires on an artefact of the correction.

    So each curve is scaled by how much variance it explains beyond the null
    expectation for fitting noise (``edf/n``). A real GC gradient keeps nearly
    all of its amplitude; a curve tracing scatter shrinks to zero.
    """
    df = counts[counts["usable"]] if usable_only else counts
    df = df[np.isfinite(df["gc"])]
    curves: dict[str, GcCurve] = {}
    for enzyme, grp in df.groupby("enzyme", sort=True):
        if len(grp) < min_anchors or grp["gc"].nunique() < 5:
            continue
        gc = grp["gc"].to_numpy(dtype=float)
        y = np.log2(grp["count"].to_numpy(dtype=float) + _PSEUDO)
        sm = lowess(y, gc, frac=frac, it=2, return_sorted=True)
        xs, ys = sm[:, 0], sm[:, 1]
        # Collapse duplicate x values so np.interp gets a strictly sorted grid.
        xs, idx = np.unique(xs, return_index=True)
        ys = ys[idx]
        if xs.size < 2:
            continue
        # Centre on the anchor-weighted mean so the correction moves anchors
        # relative to each other without shifting the enzyme's overall level
        # (which the per-enzyme intercept would absorb anyway).
        offset = ys - float(np.average(np.interp(gc, xs, ys)))

        fitted = np.interp(gc, xs, ys)
        sse = float(np.sum((y - fitted) ** 2))
        sst = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - sse / sst if sst > 0 else 0.0
        # Effective degrees of freedom of a lowess smoother, approximately
        # 1.2/frac; the R^2 it would reach on pure noise is edf/n.
        edf = min(float(len(gc)), 1.2 / max(frac, 1e-6))
        r2_null = edf / len(gc)
        shrink = float(np.clip((r2 - r2_null) / r2, 0.0, 1.0)) if r2 > 0 else 0.0
        offset = offset * shrink

        curves[str(enzyme)] = GcCurve(
            enzyme=str(enzyme),
            gc=xs,
            log2_offset=offset,
            n_anchors=int(len(grp)),
            amplitude=float(np.ptp(offset)),
            shrinkage=shrink,
            r2=float(r2),
        )
    return curves


def add_anchor_offsets(counts: pd.DataFrame, curves: dict[str, GcCurve]) -> pd.DataFrame:
    """Return ``counts`` with a ``gc_offset`` column.

    The offset is what the anchor's local GC is estimated to *add* to its log2
    count; the window rate is corrected by subtracting the window mean of it.
    Anchors of an enzyme with no fitted curve get 0.
    """
    out = counts.copy()
    off = np.zeros(len(out), dtype=float)
    for enzyme, curve in curves.items():
        sel = (out["enzyme"] == enzyme).to_numpy()
        if sel.any():
            off[sel] = curve(out.loc[sel, "gc"].to_numpy())
    out["gc_offset"] = off
    return out


def apply_to_windows(windows: pd.DataFrame, column: str = "log2_rate") -> pd.DataFrame:
    """Subtract each window's mean anchor GC offset from its log2 rate.

    Requires ``mean_gc_offset``, which :func:`sk2bgrow.ztp.window_rates` fills in
    when the count table carries ``gc_offset``. Without that column the frame is
    returned unchanged and marked ``gc_corrected = False``, so a run that skipped
    correction can never be mistaken for one that applied it.
    """
    out = windows.copy()
    if "mean_gc_offset" not in out.columns:
        out["gc_corrected"] = False
        out[f"{column}_raw"] = out[column]
        return out
    out[f"{column}_raw"] = out[column]
    out[column] = out[column] - out["mean_gc_offset"].fillna(0.0)
    out["gc_corrected"] = True
    return out


def tukey_mask(values: np.ndarray, k: float = 1.5) -> np.ndarray:
    """Boolean mask of points inside the Tukey fences ``[Q1 - k*IQR, Q3 + k*IQR]``.

    Pilea uses the same fence to drop outlier windows before fitting. NaNs are
    masked out; a degenerate (zero-IQR) distribution keeps everything, since an
    IQR of 0 would otherwise reject every value that is not the median.
    """
    v = np.asarray(values, dtype=float)
    finite = np.isfinite(v)
    if finite.sum() < 4:
        return finite
    q1, q3 = np.percentile(v[finite], [25, 75])
    iqr = q3 - q1
    if iqr <= 0:
        return finite
    return finite & (v >= q1 - k * iqr) & (v <= q3 + k * iqr)


def enzyme_efficiency(counts: pd.DataFrame, usable_only: bool = True) -> pd.DataFrame:
    """Median count per anchor for each enzyme, relative to the panel median.

    Reported for QC only. An enzyme whose efficiency is far below its peers in a
    given sample is a candidate for digestion failure or methylation sensitivity
    (report risk R1) — but, per this module's docstring, it does not by itself
    bias that enzyme's PTR.
    """
    df = counts[counts["usable"]] if usable_only else counts
    if df.empty:
        return pd.DataFrame(columns=["sample", "genome_id", "enzyme", "median_count", "mean_count", "n_anchors", "rel_efficiency"])
    g = (
        df.groupby(["sample", "genome_id", "enzyme"], sort=True)["count"]
        .agg(median_count="median", mean_count="mean", n_anchors="size")
        .reset_index()
    )
    ref = g.groupby(["sample", "genome_id"])["mean_count"].transform("median")
    g["rel_efficiency"] = np.where(ref > 0, g["mean_count"] / ref, np.nan)
    return g

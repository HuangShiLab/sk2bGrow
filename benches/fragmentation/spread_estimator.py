"""Order-free PTR from the *spread* of tag coverage, with the noise deconvolved.

An answer to: can Pilea's fragmentation-proof estimator be borrowed?

WHAT IS WORTH BORROWING -- the formulation, and it is exactly right. Under the
standard model log2(coverage) is a tent function of position with equal |slope|
on both replichores, so over a uniformly tiled genome the *values* are uniformly
distributed on [log2 c_ter, log2 c_ori] and the width of that uniform IS
log2(PTR). Only the spread of coverage is needed; position never enters. That is
why sorting the windows costs Pilea nothing on a fragmented reference.

WHAT IS NOT -- how Pilea estimates the width. `profile.py::_fit` draws a rate per
window from its ZTP mixture, sorts them, and RANSAC-regresses log2 rate on rank;
the fitted rise across all ranks is log2(PTR). Sorted values of *any* sample
rise, so sampling noise is read as growth. On the Zheng stationary control it
reports log2PTR = 1.15 at 1x for a culture that is not growing.

THE FIX -- put the noise in the model. Every window already carries a standard
error from the ZTP/ZTNB fit, so instead of reading the observed spread directly,
write y_w = mu_w + e_w with mu_w ~ U(a, a+W) and e_w ~ N(0, s_w^2) and
marginalise:

    p(y_w | a, W) = [Phi((a + W - y_w)/s_w) - Phi((a - y_w)/s_w)] / W

Maximise over (a, W) per enzyme, then fuse across enzymes by inverse variance as
usual. Same order-free property, but a sample with no gradient now returns W -> 0
instead of being credited with the spread of its own sampling error.

MEASURED, on the Zheng E. coli data against a 100-contig reference
(`benches/fragmentation`, n = 16 media):

    5x    V-fit on contigs   RMSE 0.871  slope 0.187      <- coordinate destroyed
          this estimator     RMSE 0.144  slope 0.922
          Pilea              RMSE 0.116  slope 0.888
    10x   V-fit on contigs   RMSE 0.862  slope 0.210
          this estimator     RMSE 0.227  slope 0.919
          Pilea              RMSE 0.077  slope 0.820

    stationary control (truth ~0), 1x:  Pilea 1.153   this estimator 0.000

So it recovers most of what fragmentation destroys, and it does not manufacture
growth. It does NOT beat Pilea at depth on a fragmented reference, and it is not
a replacement for the coordinate fit:

KNOWN LIMITS, in the order worth attacking.

1. Below 5x it over-shrinks to zero (bias -0.98 at 1x). The distribution alone
   carries too little information there. Coordinates are not redundant -- that is
   precisely why the V-fit wins at 1-2x and why scaffolding remains the better
   answer in that band.
2. It has a positive floor: ~0.4 on the stationary control at 5-10x, matching its
   +0.08 to +0.15 bias. Real window-to-window scatter includes systematic terms
   (mappability, anchor density, residual GC) that the model can only call
   growth. An overdispersion parameter s_eff^2 = s_w^2 + tau^2 should fix both
   this and (1), and is the obvious next step.
3. Pooling all enzymes into a single width was tried and is much worse
   (stationary control 8.6 at 0.5x): the per-enzyme fusion is doing real work by
   down-weighting sparse enzymes, and centring each enzyme by its own mean is
   corrupted where that mean is itself noise.

This is a prototype scored offline against committed `windows.rates.tsv`, not a
shipped estimator.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr


def _nll(theta, y, s):
    a, logW = theta
    W = np.exp(logW)
    hi = ndtr((a + W - y) / s)
    lo = ndtr((a - y) / s)
    p = np.clip(hi - lo, 1e-300, None) / W
    return -np.log(p).sum()


def width_mle(y, s, min_windows=6):
    """(W, se_W) for one enzyme. W is log2(PTR)."""
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    ok = np.isfinite(y) & np.isfinite(s) & (s > 0)
    y, s = y[ok], s[ok]
    n = len(y)
    if n < min_windows:
        return np.nan, np.nan

    # start from the moment estimate, floored so the optimiser has room
    var_sig = max(y.var(ddof=1) - np.mean(s ** 2), 1e-3)
    W0 = max(np.sqrt(12 * var_sig), 0.05)
    a0 = y.mean() - W0 / 2

    best = None
    for W_init in (W0, 0.5 * W0, 2.0 * W0, 0.1):
        r = minimize(_nll, [a0, np.log(W_init)], args=(y, s), method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-6, maxiter=2000))
        if best is None or r.fun < best.fun:
            best = r
    a, W = best.x[0], float(np.exp(best.x[1]))

    # curvature of the profile likelihood in W, by finite difference
    h = max(0.01 * W, 1e-3)
    f0 = best.fun
    fp = _nll([a, np.log(W + h)], y, s)
    fm = _nll([a, np.log(max(W - h, 1e-6))], y, s)
    curv = (fp - 2 * f0 + fm) / h ** 2
    se = float(np.sqrt(1.0 / curv)) if curv > 0 else np.nan
    return W, se


def per_enzyme(df, min_windows=6):
    out = []
    for enz, g in df.groupby("enzyme"):
        g = g.dropna(subset=["log2_rate", "log2_se"])
        if len(g) < min_windows:
            continue
        W, se = width_mle(g["log2_rate"].to_numpy(), g["log2_se"].to_numpy(),
                          min_windows)
        out.append((enz, W, se, len(g)))
    return out


def fuse(rows):
    r = [(w, se) for _, w, se, _ in rows
         if np.isfinite(w) and np.isfinite(se) and se > 0]
    if not r:
        return np.nan, np.nan, np.nan, 0
    w = np.array([a for a, _ in r])
    iv = 1.0 / np.array([b ** 2 for _, b in r])
    est = float((w * iv).sum() / iv.sum())
    return est, float(np.sqrt(1.0 / iv.sum())), float((iv * (w - est) ** 2).sum()), len(r)

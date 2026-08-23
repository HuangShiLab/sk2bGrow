"""Window rate models: zero-truncated Poisson mixtures and a negative-binomial
branch.

Reference: Pilea's ``ztp.py``. Pilea fits each window's k-mer counts to a
zero-truncated Poisson mixture (EM for parameters, BIC for the component count)
and takes the highest-weight component's rate as the window's expected coverage.
Truncation is the point: systematic sequence divergence makes some reference
k-mers genuinely absent from a sample, and counting those as true zeros drags a
plain Poisson fit downward.

Two changes the design report asks for (§7.1 step 3):

* **a negative-binomial branch.** Pilea's own dispersion filter shows residual
  overdispersion is real; NB accommodates it instead of discarding the sample.
  At low coverage NB is unidentifiable, so the model falls back to ZTP — chosen
  by BIC rather than by a hard coverage threshold.
* **per-enzyme stratification.** Windows are built inside one enzyme's anchor
  series, not across the union, so each enzyme is an independent measurement
  channel with its own efficiency and GC context.

Every rate is returned with a standard error. That is what makes the
inverse-variance fusion in :mod:`sk2bgrow.fusion` possible, and it is the piece
Pilea's bootstrap has to approximate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize, special

__all__ = [
    "MIN_ALPHA",
    "ztp_mean",
    "ztp_var",
    "solve_ztp_lambda",
    "auto_window_size",
    "ZtpMixture",
    "fit_ztp_mixture",
    "ZtnbFit",
    "fit_ztnb",
    "WindowRate",
    "estimate_window_rate",
    "window_rates",
]

_EPS = 1e-12


# --------------------------------------------------------------------------
# zero-truncated Poisson primitives
# --------------------------------------------------------------------------

def ztp_logpmf(k: np.ndarray, lam: float) -> np.ndarray:
    """log P(X = k | X >= 1) for a Poisson with rate ``lam``."""
    lam = max(float(lam), _EPS)
    k = np.asarray(k, dtype=float)
    return k * np.log(lam) - lam - special.gammaln(k + 1.0) - np.log1p(-np.exp(-lam))


def ztp_mean(lam: float) -> float:
    """E[X | X >= 1] = lam / (1 - exp(-lam)); tends to 1 as lam -> 0."""
    lam = float(lam)
    if lam < 1e-8:
        return 1.0 + lam / 2.0
    return lam / (1.0 - np.exp(-lam))


def ztp_var(lam: float) -> float:
    """Var[X | X >= 1] = m (1 + lam - m)."""
    m = ztp_mean(lam)
    return max(m * (1.0 + float(lam) - m), _EPS)


def ztp_dmean_dlam(lam: float) -> float:
    """d/dlam of :func:`ztp_mean`, for the delta-method standard error."""
    lam = float(lam)
    if lam < 1e-8:
        return 0.5
    u = 1.0 - np.exp(-lam)
    return (u - lam * np.exp(-lam)) / (u * u)


def solve_ztp_lambda(observed_mean: float, upper: float = 1e6, tol: float = 1e-10) -> float:
    """Invert :func:`ztp_mean`.

    ``ztp_mean`` is strictly increasing from 1 (at lam -> 0) upward, so a
    truncated sample mean at or below 1 carries no information about lam beyond
    "very small" and maps to ~0. Returning 0 there rather than raising keeps a
    single all-ones window from aborting a whole genome.

    Solved by safeguarded Newton using the analytic derivative
    :func:`ztp_dmean_dlam`. This sits in the innermost EM loop and is called
    ~10^5 times for one sample; Newton converges in a handful of steps against
    the dozens of function evaluations a bracketing solver needs. Brentq remains
    the fallback for the rare case Newton escapes its bracket.
    """
    m = float(observed_mean)
    if not np.isfinite(m) or m <= 1.0 + 1e-9:
        return 0.0
    # For a large truncated mean the truncation is negligible and lam ~= m; near
    # the lower limit m -> 1 + lam/2, so lam ~= 2(m-1).
    lam = m if m > 2.0 else 2.0 * (m - 1.0)
    lam = min(max(lam, 1e-9), upper)
    for _ in range(40):
        f = ztp_mean(lam) - m
        if abs(f) < tol * max(1.0, m):
            return float(lam)
        d = ztp_dmean_dlam(lam)
        if d <= 0 or not np.isfinite(d):
            break
        step = f / d
        new_lam = lam - step
        if not np.isfinite(new_lam) or new_lam <= 0:
            new_lam = lam / 2.0
        if new_lam > upper:
            return upper
        if abs(new_lam - lam) < tol * max(1.0, lam):
            return float(new_lam)
        lam = new_lam
    # Fallback: bracket and bisect.
    lo, hi = 1e-9, max(m, 1.0)
    while ztp_mean(hi) < m:
        hi *= 2.0
        if hi > upper:
            return upper
    return float(optimize.brentq(lambda x: ztp_mean(x) - m, lo, hi, xtol=1e-10, rtol=1e-12))


# --------------------------------------------------------------------------
# ZTP mixture, fitted by EM with BIC model selection
# --------------------------------------------------------------------------

@dataclass
class ZtpMixture:
    """A fitted mixture of zero-truncated Poissons."""

    lambdas: np.ndarray
    weights: np.ndarray
    loglik: float
    bic: float
    n_obs: int
    n_iter: int
    converged: bool
    #: Responsibility mass of each component, i.e. its effective sample size.
    n_eff: np.ndarray = field(default_factory=lambda: np.array([]))

    @property
    def n_components(self) -> int:
        return int(len(self.lambdas))

    @property
    def dominant(self) -> int:
        """Index of the highest-weight component — Pilea's choice of rate."""
        return int(np.argmax(self.weights))

    @property
    def rate(self) -> float:
        return float(self.lambdas[self.dominant])


def _mixture_loglik(counts: np.ndarray, lambdas: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray]:
    comp = np.stack([np.log(max(w, _EPS)) + ztp_logpmf(counts, lam) for w, lam in zip(weights, lambdas)])
    total = special.logsumexp(comp, axis=0)
    resp = np.exp(comp - total)
    return float(total.sum()), resp


def fit_ztp_mixture(
    counts: np.ndarray,
    n_components: int = 1,
    max_iter: int = 200,
    tol: float = 1e-7,
    seed: int = 0,
) -> ZtpMixture:
    """EM fit of a ``n_components`` ZTP mixture to positive counts.

    ``counts`` must already be truncated: zeros are dropped by the caller, since
    "which anchors were never seen" is a different question (containment) from
    "how deep are the ones that were".
    """
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts >= 1]
    n = counts.size
    if n == 0:
        return ZtpMixture(np.array([0.0]), np.array([1.0]), -np.inf, np.inf, 0, 0, False, np.array([0.0]))

    if n_components == 1:
        # Single component: the MLE is one inversion of the mean, no EM at all.
        lam = solve_ztp_lambda(float(counts.mean()))
        lambdas = np.array([lam])
        weights = np.array([1.0])
        loglik = float(ztp_logpmf(counts, lam).sum()) if lam > 0 else -np.inf
        return ZtpMixture(
            lambdas=lambdas,
            weights=weights,
            loglik=loglik,
            bic=float(-2.0 * loglik + np.log(n)),
            n_obs=n,
            n_iter=1,
            converged=True,
            n_eff=np.array([float(n)]),
        )

    # Initialise components at quantiles of the data, which spreads them along
    # the observed range.
    rng = np.random.default_rng(seed)
    qs = np.linspace(0.15, 0.85, n_components)
    starts = np.quantile(counts, qs) if n > 1 else np.repeat(counts[0], n_components)
    lambdas = np.array([solve_ztp_lambda(max(s, 1.0 + 1e-6)) for s in starts], dtype=float)
    lambdas = np.where(lambdas <= 0, rng.uniform(0.1, 1.0, size=n_components), lambdas)
    lambdas.sort()
    weights = np.full(n_components, 1.0 / n_components)

    prev = -np.inf
    resp = np.ones((n_components, n)) / n_components
    it = 0
    converged = False
    for it in range(1, max_iter + 1):
        loglik, resp = _mixture_loglik(counts, lambdas, weights)
        if np.isfinite(prev) and abs(loglik - prev) < tol * max(1.0, abs(prev)):
            converged = True
            break
        prev = loglik
        mass = resp.sum(axis=1)
        weights = np.maximum(mass / n, _EPS)
        weights /= weights.sum()
        for j in range(n_components):
            if mass[j] <= _EPS:
                continue
            lambdas[j] = solve_ztp_lambda((resp[j] * counts).sum() / mass[j])

    loglik, resp = _mixture_loglik(counts, lambdas, weights)
    n_params = 2 * n_components - 1
    bic = -2.0 * loglik + n_params * np.log(n)
    return ZtpMixture(
        lambdas=lambdas,
        weights=weights,
        loglik=loglik,
        bic=float(bic),
        n_obs=n,
        n_iter=it,
        converged=converged,
        n_eff=resp.sum(axis=1),
    )


# --------------------------------------------------------------------------
# zero-truncated negative binomial
# --------------------------------------------------------------------------

@dataclass
class ZtnbFit:
    """A zero-truncated negative binomial, parameterised by mean and dispersion."""

    mu: float
    #: Overdispersion: Var = mu + alpha * mu^2. alpha -> 0 recovers Poisson.
    alpha: float
    loglik: float
    bic: float
    n_obs: int

    @property
    def rate(self) -> float:
        """The untruncated mean, comparable to a ZTP ``lambda``."""
        return self.mu


#: Smallest dispersion the NB branch may report. Below this the model is
#: numerically indistinguishable from Poisson, and letting the optimiser chase
#: alpha -> 0 only finds rounding noise.
MIN_ALPHA = 1e-8


def _log_rising(k: np.ndarray, r: float) -> np.ndarray:
    """``lgamma(k + r) - lgamma(r)`` for non-negative integer ``k``.

    Written as the cumulative sum ``sum_i log(r + i)`` rather than a difference
    of two ``gammaln`` calls. When alpha is small, ``r = 1/alpha`` is huge and
    both gammaln terms are ~1e13 while their difference is ~10 — the subtraction
    then loses most of the mantissa and hands the optimiser a dozen units of
    phantom likelihood, which makes NB "win" on perfectly Poisson data.
    """
    kmax = int(k.max()) if k.size else 0
    if kmax > 1_000_000:  # pathological depth: fall back rather than allocate
        return special.gammaln(k + r) - special.gammaln(r)
    table = np.concatenate([[0.0], np.cumsum(np.log(r + np.arange(kmax)))])
    return table[k.astype(np.int64)]


def _ztnb_logpmf(k: np.ndarray, mu: float, alpha: float) -> np.ndarray:
    mu = max(float(mu), _EPS)
    alpha = max(float(alpha), MIN_ALPHA)
    r = 1.0 / alpha
    # log p and log(1-p) via log1p, so neither cancels for large r.
    log_p = -np.log1p(mu / r)
    log_q = np.log(mu) - np.log(r) - np.log1p(mu / r)
    log_nb = _log_rising(np.asarray(k), r) - special.gammaln(k + 1.0) + r * log_p + k * log_q
    log_zero = r * log_p                        # log P(X = 0)
    return log_nb - np.log1p(-np.exp(min(log_zero, -_EPS)))


def fit_ztnb(counts: np.ndarray) -> ZtnbFit:
    """MLE fit of a zero-truncated negative binomial to positive counts."""
    counts = np.asarray(counts, dtype=float)
    counts = counts[counts >= 1]
    n = counts.size
    if n == 0:
        return ZtnbFit(0.0, 0.0, -np.inf, np.inf, 0)

    def nll(theta: np.ndarray) -> float:
        mu, alpha = np.exp(theta)
        if not np.isfinite(mu) or not np.isfinite(alpha):
            return 1e18
        val = -_ztnb_logpmf(counts, mu, alpha).sum()
        return float(val) if np.isfinite(val) else 1e18

    m = float(counts.mean())
    v = float(counts.var(ddof=1)) if n > 1 else m
    alpha0 = max((v - m) / max(m * m, _EPS), 1e-3)
    x0 = np.log([max(m, 1e-3), min(max(alpha0, 1e-3), 10.0)])
    res = optimize.minimize(
        nll,
        x0,
        method="L-BFGS-B",
        bounds=[(np.log(1e-6), np.log(1e6)), (np.log(MIN_ALPHA), np.log(1e3))],
    )
    mu, alpha = np.exp(res.x)
    loglik = -float(res.fun)
    return ZtnbFit(mu=float(mu), alpha=float(alpha), loglik=loglik, bic=float(-2.0 * loglik + 2.0 * np.log(n)), n_obs=n)


def ztnb_mean(mu: float, alpha: float) -> float:
    """E[X | X >= 1] for a zero-truncated negative binomial."""
    r = 1.0 / max(alpha, MIN_ALPHA)
    p0 = np.exp(-r * np.log1p(mu / r))
    return float(mu / max(1.0 - p0, _EPS))


def ztnb_var(mu: float, alpha: float) -> float:
    """Var[X | X >= 1] for a zero-truncated negative binomial."""
    r = 1.0 / max(alpha, MIN_ALPHA)
    p0 = np.exp(-r * np.log1p(mu / r))
    s = max(1.0 - p0, _EPS)
    untrunc_var = mu + alpha * mu * mu
    m = mu / s
    # Var = (Var_untrunc + mu^2)/s - m^2, i.e. E[X^2|X>=1] - E[X|X>=1]^2.
    return float(max((untrunc_var + mu * mu) / s - m * m, _EPS))


# --------------------------------------------------------------------------
# the public per-window estimate
# --------------------------------------------------------------------------

@dataclass
class WindowRate:
    """Expected per-anchor count in one window, with its uncertainty."""

    rate: float
    se: float
    model: str
    #: Index of dispersion (Var/mean) of the positive counts; > 1 means the
    #: window is overdispersed relative to Poisson.
    dispersion: float
    n_anchors: int
    n_positive: int
    #: Fraction of the window's anchors that were observed at all.
    detected_fraction: float
    n_components: int = 1
    bic: float = np.nan

    @property
    def log2_rate(self) -> float:
        return float(np.log2(self.rate)) if self.rate > 0 else np.nan

    @property
    def log2_se(self) -> float:
        """Delta-method standard error on the log2 scale."""
        if self.rate <= 0 or not np.isfinite(self.se):
            return np.nan
        return float(self.se / (self.rate * np.log(2.0)))


def estimate_window_rate(
    counts: np.ndarray,
    model: str = "auto",
    max_components: int = 3,
) -> WindowRate:
    """Estimate one window's expected per-anchor count.

    ``model`` is ``"ztp"``, ``"nb"`` or ``"auto"``. ``"auto"`` fits both and
    keeps the lower BIC, which lets a shallow window fall back to ZTP without a
    hand-set coverage threshold.

    The standard error comes from the delta method::

        se(rate) = sqrt(Var_model[X | X >= 1] / n_eff) / |d mean / d rate|

    i.e. the sampling error of the truncated mean, propagated through the
    truncation. It is the same construction for every model, so BIC can switch
    between them without making the errors incomparable.
    """
    counts = np.asarray(counts, dtype=float)
    n_anchors = int(counts.size)
    pos = counts[counts >= 1]
    n_pos = int(pos.size)
    detected = n_pos / n_anchors if n_anchors else np.nan
    disp = float(pos.var(ddof=1) / pos.mean()) if n_pos > 1 and pos.mean() > 0 else np.nan

    if n_pos == 0:
        return WindowRate(np.nan, np.nan, "empty", disp, n_anchors, 0, detected)

    best_ztp: ZtpMixture | None = None
    if model in ("ztp", "auto"):
        # More components than distinct values cannot be identified.
        cap = max(1, min(max_components, n_pos // 10, int(np.unique(pos).size)))
        for j in range(1, cap + 1):
            fit = fit_ztp_mixture(pos, n_components=j)
            if best_ztp is None or fit.bic < best_ztp.bic:
                best_ztp = fit

    best_nb: ZtnbFit | None = None
    # NB needs both enough observations and visible overdispersion to be
    # identifiable; forcing it on a handful of counts yields a huge alpha and a
    # meaningless rate.
    if model == "nb" or (model == "auto" and n_pos >= 10):
        best_nb = fit_ztnb(pos)

    use_nb = best_nb is not None and (best_ztp is None or best_nb.bic < best_ztp.bic)

    if use_nb:
        assert best_nb is not None
        mu, alpha = best_nb.mu, best_nb.alpha
        var = ztnb_var(mu, alpha)
        h = max(mu * 1e-5, 1e-8)
        dmean = (ztnb_mean(mu + h, alpha) - ztnb_mean(max(mu - h, _EPS), alpha)) / (2.0 * h)
        se = np.sqrt(var / n_pos) / max(abs(dmean), _EPS)
        return WindowRate(
            rate=float(mu),
            se=float(se),
            model="ztnb",
            dispersion=disp,
            n_anchors=n_anchors,
            n_positive=n_pos,
            detected_fraction=detected,
            n_components=1,
            bic=best_nb.bic,
        )

    assert best_ztp is not None
    j = best_ztp.dominant
    lam = float(best_ztp.lambdas[j])
    n_eff = float(best_ztp.n_eff[j]) if best_ztp.n_eff.size else float(n_pos)
    if lam <= 0 or n_eff <= 0:
        return WindowRate(np.nan, np.nan, "ztp", disp, n_anchors, n_pos, detected, best_ztp.n_components, best_ztp.bic)
    se = np.sqrt(ztp_var(lam) / n_eff) / max(abs(ztp_dmean_dlam(lam)), _EPS)
    return WindowRate(
        rate=lam,
        se=float(se),
        model="ztp",
        dispersion=disp,
        n_anchors=n_anchors,
        n_positive=n_pos,
        detected_fraction=detected,
        n_components=best_ztp.n_components,
        bic=best_ztp.bic,
    )


#: Windows to aim for per enzyme in ``anchors_per_window="auto"`` mode. A line
#: through ~25 points is well conditioned; far fewer and the slope rides on a
#: handful of windows.
TARGET_WINDOWS_PER_ENZYME = 25

#: Never cut a window below this many anchors — a rate estimate needs a sample.
MIN_ANCHORS_PER_WINDOW = 25


def auto_window_size(n_anchors: int, cap: int = 100) -> int:
    """Anchors per window for one enzyme, given how many anchors it has.

    A single fixed size cannot serve the panel. Enzyme densities span 20x
    (report §4.1: CjeI 1 910/Mb, PpiI 73/Mb), so at a flat 100 anchors/window on
    E. coli, CjeI gets 88 windows while PpiI, BplI and PsrI get 3-4 — below any
    sensible minimum, so the three sparsest enzymes drop out of a 16-enzyme
    design entirely. Sizing per enzyme keeps all sixteen channels alive: sparse
    enzymes get smaller, noisier windows and *earn less weight* in the
    inverse-variance fusion, which is the right way for them to count less.
    """
    return int(min(max(n_anchors // TARGET_WINDOWS_PER_ENZYME, MIN_ANCHORS_PER_WINDOW), cap))


def window_rates(
    counts: pd.DataFrame,
    anchors_per_window: int | str = "auto",
    model: str = "auto",
    usable_only: bool = True,
    use_precomputed_windows: bool = False,
    window_cap: int = 100,
) -> pd.DataFrame:
    """Window a count table and estimate a rate for every window.

    Windows are cut **inside each (genome, enzyme) series**, in position order,
    with a fixed anchor count. That is the report's §6.3 adaptive window: equal
    statistical power per window, rather than equal base pairs. Cutting per
    enzyme rather than across the 16-enzyme union is what keeps each enzyme an
    independent channel for :mod:`sk2bgrow.fusion`.

    ``anchors_per_window`` is an integer, or ``"auto"`` to size each enzyme by
    :func:`auto_window_size` (capped at ``window_cap``).

    Set ``use_precomputed_windows`` to group by the Rust ``window_id`` instead —
    that reproduces Pilea's fixed-bp windowing for an A/B benchmark.

    Windows never span contigs: two anchors on different contigs have no defined
    distance, so pooling them would fabricate a coordinate.
    """
    df = counts
    if usable_only:
        df = df[df["usable"]]
    if df.empty:
        return pd.DataFrame(
            columns=[
                "sample", "genome_id", "genome", "enzyme", "window", "contig_id",
                "start", "end", "global_mid", "n_anchors", "n_positive",
                "anchors_per_window", "detected_fraction", "rate", "se", "log2_rate", "log2_se",
                "model", "dispersion", "mean_gc", "mean_gc_offset",
            ]
        )

    rows = []
    group_cols = ["sample", "genome_id", "genome", "enzyme"]
    for keys, grp in df.groupby(group_cols, sort=True):
        grp = grp.sort_values(["contig_id", "position"])
        n_per_window = (
            auto_window_size(len(grp), window_cap)
            if isinstance(anchors_per_window, str)
            else int(anchors_per_window)
        )
        if use_precomputed_windows:
            grp = grp[grp["window_id"] != 0xFFFFFFFF]
            blocks = list(grp.groupby("window_id", sort=True))
        else:
            blocks = []
            for contig, cg in grp.groupby("contig_id", sort=True):
                n = len(cg)
                if n == 0:
                    continue
                # A trailing block far below target is folded into the previous
                # window rather than reported with inflated variance.
                idx = np.arange(n) // n_per_window
                if n % n_per_window and n // n_per_window >= 1:
                    tail = n % n_per_window
                    if tail * 2 < n_per_window:
                        idx[idx == idx.max()] = idx.max() - 1
                for w, wg in cg.groupby(idx, sort=True):
                    blocks.append(((contig, int(w)), wg))

        for wkey, wg in blocks:
            est = estimate_window_rate(wg["count"].to_numpy(), model=model)
            sample, genome_id, genome, enzyme = keys
            rows.append(
                {
                    "sample": sample,
                    "genome_id": int(genome_id),
                    "genome": genome,
                    "enzyme": enzyme,
                    "window": f"{wkey}",
                    "contig_id": int(wg["contig_id"].iloc[0]),
                    "start": int(wg["position"].min()),
                    "end": int(wg["position"].max()),
                    "global_mid": float(wg["global_position"].mean()),
                    "anchors_per_window": n_per_window,
                    "n_anchors": est.n_anchors,
                    "n_positive": est.n_positive,
                    "detected_fraction": est.detected_fraction,
                    "rate": est.rate,
                    "se": est.se,
                    "log2_rate": est.log2_rate,
                    "log2_se": est.log2_se,
                    "model": est.model,
                    "dispersion": est.dispersion,
                    "mean_gc": float(np.nanmean(wg["gc"].to_numpy())) if "gc" in wg else np.nan,
                    # Mean anchor-level GC offset, so gc_bias.apply_to_windows can
                    # correct at anchor resolution rather than on the window mean GC.
                    "mean_gc_offset": (
                        float(np.nanmean(wg["gc_offset"].to_numpy())) if "gc_offset" in wg else np.nan
                    ),
                }
            )
    return pd.DataFrame(rows)

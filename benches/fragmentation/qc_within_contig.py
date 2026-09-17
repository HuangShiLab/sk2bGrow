#!/usr/bin/env python3
"""R3: a QC statistic that can see a destroyed coordinate.

THE IDEA. The current QC asks whether the sixteen enzyme strata agree
(Cochran's Q). A scrambled reference makes all sixteen agree there is no
gradient, so Q passes 100% of fragmented estimates. But scrambling only
destroys the coordinate *between* contigs. Within a contig the sequence --
and therefore the local replication gradient -- is untouched. Scrambling
also breaks continuity *at* contig boundaries: window rates stop matching
across a join. Two statistics exploit those two survivors:

  b_hat   the pipeline's fused V-fit amplitude (uses the global coordinate;
          destroyed by scrambling),

  A_loc   a within-contig amplitude: the RMS of per-contig slopes, scaled by
          half the genome length, with the sampling-noise floor subtracted
          (uses only coordinates *inside* each contig; indifferent to order
          and orientation, so scrambling cannot hurt it). It needs several
          windows per contig per enzyme to separate slope from enzyme
          offsets, so it is trusted only where contigs are long enough
          (this panel: N up to ~20; see R3_QC.md).

  J       the boundary-jump excess: mean squared step between consecutive
          windows across a contig boundary, minus the same quantity within
          contigs, both corrected for the window SEs. A correct coordinate
          is continuous, so J ~ 0; a scrambled one jumps to an unrelated
          coverage level at every boundary, so J ~ log2PTR^2 / 6. It needs
          many boundaries, so it strengthens exactly where A_loc weakens.

The four regimes separate:

  complete / correctly ordered, growing:   A_loc ~ b_hat, J ~ 0    -> quiet
  scrambled, growing:                      A_loc >> b_hat or J >> 0 -> FIRE
  stationary (log2PTR = 0), any reference: A_loc ~ 0, J ~ 0        -> quiet
  complete (N = 1 contig):                 both undefined          -> quiet

Neither statistic can fire on a stationary culture: there is no local
gradient for the global fit to lose, and there are no boundary jumps to
find. That is exactly the distinction Cochran's Q cannot make.

DEFINITIONS. Enzymes are median-centred genome-wide (as the origin search
already does). For each contig c, weighted least squares over its windows
(weight 1/log2_se^2), x in bp contig-local:

    y_w = a_c + m_c * (x_w - xbar_c)

    A_loc = sqrt(max(0, sum w_c m_c^2 / sum w_c  -  k / sum w_c)) * (L / 2)

with w_c = 1/s_c^2, k the number of usable contigs, L the total reference
length. k/sum(w_c) is E[Q] under "no slopes", so A_loc is noise-corrected:
a stationary sample returns ~0 rather than the spread of its own noise.

    z_loc   = A_loc / se(A_loc)                      is there a local gradient?
    z_short = (A_loc - b_hat) / sqrt(se_A^2 + se_b^2)  did the global fit lose it?

For J, within each enzyme series ordered by reference coordinate, take
consecutive window pairs with step d = y2 - y1 and noise floor
f = s1^2 + s2^2; the excess e = d^2 - f estimates the squared true step.
J is the mean excess of boundary-crossing pairs minus the mean excess of
within-contig pairs, with an empirical (per-pair scatter) standard error,
and z_J = J / se(J).

The slope branch is used only where contigs carry enough within-contig
information to separate slope from enzyme offset: median windows per usable
contig >= SLOPE_MIN_WINDOWS. The jump branch requires >= JUMP_MIN_PAIRS
boundary pairs. The proposed gate fires when

    slope branch:  usable AND z_loc > 2 AND z_short > 2 AND A_loc - b_hat > 0.5
    jump branch:   usable AND z_J > 3.5 AND J > 0.03

(the magnitude margins matter: at 1x the V-fit itself runs ~25% low -- the
A4 compression -- and a pure significance test cannot tell that from a
destroyed coordinate).

Prototype scored offline against the windows.rates.tsv / output.tsv that
profile runs already write; not wired into the pipeline.
"""
import numpy as np
import pandas as pd

MIN_WINDOWS_PER_CONTIG = 8
SLOPE_MIN_WINDOWS = 20   # median per usable contig; ~1.25 windows/enzyme
JUMP_MIN_PAIRS = 8


def per_contig_slopes(windows, min_windows=MIN_WINDOWS_PER_CONTIG):
    """One slope (log2 rate per bp) with SE per contig.

    Enzymes are first median-centred genome-wide (as the pipeline's origin
    search does), so each contig then needs only a 2-parameter fit
    (intercept + slope) instead of one intercept per enzyme. That matters:
    at N ~ 100 contigs a contig holds only ~2 windows per enzyme, and a
    per-enzyme-intercept fit is near-saturated there -- its slopes are
    unstable and their standard errors are not honest. The contig intercept
    stays free: a scrambled reference moves contigs to arbitrary coverage
    levels.
    """
    df = windows.dropna(subset=["log2_rate", "log2_se"])
    df = df[df["log2_se"] > 0]
    if df.empty:
        return pd.DataFrame()
    centre = df.groupby("enzyme")["log2_rate"].median()
    df = df.assign(y=df["log2_rate"] - df["enzyme"].map(centre))

    out = []
    for cid, g in df.groupby("contig_id"):
        if len(g) < min_windows or g["enzyme"].nunique() < 2:
            continue
        x = 0.5 * (g["start"].to_numpy(float) + g["end"].to_numpy(float))
        if x.max() - x.min() < 1e4:  # no within-contig coordinate to fit
            continue
        y = g["y"].to_numpy(float)
        w = 1.0 / g["log2_se"].to_numpy(float) ** 2
        xc = x - np.average(x, weights=w)
        X = np.column_stack([np.ones_like(xc), xc])
        Xw = X * w[:, None]
        with np.errstate(all="ignore"):
            XtX_inv = np.linalg.pinv(Xw.T @ X)
            beta = XtX_inv @ (Xw.T @ y)
        resid = y - X @ beta
        chi2 = float((w * resid ** 2).sum())
        dof = len(g) - 2
        # The window SEs are not always calibrated (the ZTP/ZTNB layer at
        # near-zero counts is the A4 suspect). Scale the slope SE by the
        # contig's own residual chi2 when it exceeds its degrees of freedom,
        # so a contig whose windows scatter more than their SEs claim does
        # not get an artificially tight slope.
        scale = max(1.0, chi2 / dof) if dof > 0 else np.inf
        se = float(np.sqrt(max(XtX_inv[1, 1], 1e-30) * scale))
        if not (np.isfinite(beta[1]) and np.isfinite(se)):
            continue
        out.append(dict(contig_id=int(cid), slope=beta[1],
                        slope_se=se,
                        n_windows=len(g), n_enzymes=g["enzyme"].nunique(),
                        span=float(x.max() - x.min())))
    return pd.DataFrame(out)


def local_amplitude(slopes, total_len):
    """Noise-corrected RMS within-contig slope, rescaled to log2PTR units.

    Returns (A_loc, se_A, k usable contigs). Undefined (NaN) for k < 2.
    """
    if len(slopes) < 2:
        return np.nan, np.nan, len(slopes)
    m = slopes["slope"].to_numpy(float)
    s = slopes["slope_se"].to_numpy(float)
    ok = np.isfinite(m) & np.isfinite(s) & (s > 0)
    m, s = m[ok], s[ok]
    k = len(m)
    if k < 2:
        return np.nan, np.nan, k
    w = 1.0 / s ** 2
    W = w.sum()
    Q = (w * m ** 2).sum() / W
    floor = k / W
    # Var(m^2) = 2 s^4 + 4 m^2 s^2 for m ~ N(mu, s^2); propagate to Q.
    var_Q = ((2.0 + 4.0 * m ** 2 / s ** 2) / W ** 2).sum()
    se_Q = float(np.sqrt(var_Q))
    A2 = max(Q - floor, 0.0)
    A = float(np.sqrt(A2))
    se_A = se_Q / (2.0 * A) if A > 0 else float(np.sqrt(max(se_Q, 1e-30)) / 2)
    return A * (total_len / 2.0), se_A * (total_len / 2.0), k


def boundary_jumps(windows, min_pairs=JUMP_MIN_PAIRS):
    """Boundary-jump excess J: do window rates leap at contig boundaries?

    Within each enzyme series (ordered by reference coordinate), consecutive
    window pairs give a step d = y2 - y1 with noise floor f = s1^2 + s2^2.
    The excess e = d^2 - f estimates the squared *true* step. J is the mean
    excess of boundary-crossing pairs minus within-contig pairs: ~0 on any
    continuous coordinate (correct assembly, or no gradient at all), about
    log2PTR^2/6 on a scrambled one.

    Returns (J, se_J, n_boundary_pairs). SE is empirical from the per-pair
    scatter, so it does not trust the window SEs beyond the floor subtraction.
    """
    df = windows.dropna(subset=["log2_rate", "log2_se"])
    df = df[df["log2_se"] > 0]
    # A window with log2_se > 1 constrains its rate to less than a factor of
    # two; pairs involving it carry no usable step information, only noise.
    df = df[df["log2_se"] <= 1.0]
    e_b, e_w = [], []
    for _, g in df.groupby("enzyme"):
        g = g.sort_values(["contig_id", "start"])
        y = g["log2_rate"].to_numpy(float)
        f = g["log2_se"].to_numpy(float) ** 2
        same = np.diff(g["contig_id"].to_numpy()) == 0
        d2 = np.diff(y) ** 2
        f = f[:-1] + f[1:]
        # Winsorize the squared step at 9x its noise floor (3 sigma): one
        # wild low-count window must not dominate a mean.
        e = np.minimum(d2, 9.0 * f) - f
        e_b.append(e[~same])
        e_w.append(e[same])
    e_b = np.concatenate(e_b) if e_b else np.array([])
    e_w = np.concatenate(e_w) if e_w else np.array([])
    if len(e_b) < min_pairs or len(e_w) < min_pairs:
        return np.nan, np.nan, len(e_b)
    J = float(e_b.mean() - e_w.mean())
    se = float(np.sqrt(e_b.var(ddof=1) / len(e_b) + e_w.var(ddof=1) / len(e_w)))
    return J, se, len(e_b)


def wcg_stat(windows, total_len, b_hat, se_b,
             slope_min_windows=SLOPE_MIN_WINDOWS):
    """The within-contig gradient (WCG) gate for one profile run.

    Returns both branches (slope-amplitude shortfall, boundary jumps) and
    the proposed combined flag.
    """
    slopes = per_contig_slopes(windows)
    A, se_A, k = local_amplitude(slopes, total_len)
    J, se_J, n_b = boundary_jumps(windows)
    out = dict(k_contigs=k, A_loc=A, se_A=se_A, b_hat=b_hat, se_b=se_b,
               slope_usable=False, z_loc=np.nan, z_short=np.nan,
               J=J, se_J=se_J, n_boundary=n_b, z_J=np.nan,
               flag_slope=False, flag_jump=False, flag=False)

    if np.isfinite(A) and se_A > 0 and len(slopes) >= 2:
        out["slope_usable"] = bool(
            slopes["n_windows"].median() >= slope_min_windows)
        out["z_loc"] = A / se_A
        if np.isfinite(b_hat) and np.isfinite(se_b) and se_b > 0:
            out["z_short"] = (A - b_hat) / float(np.sqrt(se_A ** 2 + se_b ** 2))
        out["flag_slope"] = bool(
            out["slope_usable"] and out["z_loc"] > 2.0
            and np.isfinite(out["z_short"]) and out["z_short"] > 2.0
            and A - b_hat > 0.5)
    if np.isfinite(J) and se_J > 0:
        out["z_J"] = J / se_J
        out["flag_jump"] = bool(out["z_J"] > 3.5 and J > 0.03)
    out["flag"] = bool(out["flag_slope"] or out["flag_jump"])
    return out

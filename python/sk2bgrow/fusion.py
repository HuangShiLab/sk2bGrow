"""Cross-enzyme fusion and consistency QC — the project's core new step.

Each of the 16 Type IIB enzymes measures the *same* biological quantity (the
ori-ter copy-number gradient) through a *different* set of loci, with its own
digestion efficiency and its own GC neighbourhood. That is a stratified design
with real replication built in, and it is the thing no single random sketch can
offer (report defect D4: Pilea's uncertainty comes from resampling mixture
components, which is a computational pseudo-replicate, not an independent
measurement).

Two products come out of that structure:

**A fused estimate.** Inverse-variance weighting is the minimum-variance linear
combination of unbiased estimates, so enzymes with more anchors and cleaner fits
count for more, automatically.

**A consistency test.** Under the null that every enzyme measures one common
value, Cochran's Q is chi-square distributed on ``k-1`` degrees of freedom.
A significant Q means the enzymes *disagree* — an enzyme with too few anchors, a
methylation-blocked site class, a mis-assembled region. That is a genuine QC
signal available at zero extra sequencing cost.

When Q rejects, the fixed-effect standard error is known to be too small: it
assumes the only scatter is sampling noise. The estimator then escalates to the
DerSimonian-Laird random-effects weights, which add the between-enzyme variance
component. Reporting a tight interval around a value the enzymes visibly
disagree about would be the worst of both worlds.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

__all__ = ["FusionResult", "fuse", "fuse_table", "MIN_ANCHORS_PER_ENZYME"]

#: Below this many anchors an enzyme's estimate is treated as uninformative and
#: excluded before fusion. The report's density table shows PpiI and PsrI at
#: 70-140 anchors/Mb, so on a small genome they legitimately fall under it.
MIN_ANCHORS_PER_ENZYME = 30

#: Q-test p-value below which the enzymes are considered inconsistent.
DEFAULT_ALPHA = 0.05


@dataclass
class FusionResult:
    """A fused PTR estimate plus the evidence for trusting it."""

    log2_ptr: float
    se: float
    n_enzymes: int
    #: Cochran's Q and its p-value under the homogeneity null.
    q: float = np.nan
    q_pvalue: float = np.nan
    #: I^2: the share of total variation attributable to real between-enzyme
    #: differences rather than sampling error.
    i2: float = np.nan
    #: DerSimonian-Laird between-enzyme variance.
    tau2: float = 0.0
    model: str = "fixed"
    #: Enzymes dropped, with the reason.
    excluded: dict[str, str] = field(default_factory=dict)
    used: tuple[str, ...] = ()
    #: Per-enzyme standardised residuals, for spotting the culprit.
    residuals: dict[str, float] = field(default_factory=dict)
    ok: bool = True
    note: str = ""

    @property
    def ptr(self) -> float:
        return float(2.0**self.log2_ptr) if np.isfinite(self.log2_ptr) else np.nan

    @property
    def consistency_checked(self) -> bool:
        """True when Q could actually be computed, i.e. at least two enzymes.

        With one surviving enzyme there are zero degrees of freedom, so
        ``q_pvalue`` is NaN and :attr:`consistent` reports ``True`` for a test
        that never ran. Filter on ``consistency_checked and consistent`` when
        what you mean is "the enzymes were checked and agreed".
        """
        return self.n_enzymes >= 2 and bool(np.isfinite(self.q_pvalue))

    @property
    def consistent(self) -> bool:
        """True when the enzymes do not measurably disagree.

        **Also True when no check was possible** (a single surviving enzyme
        gives NaN ``q_pvalue``). That is deliberate — one enzyme's estimate is
        not wrong merely for being alone — but it makes this property unsafe as
        a QC filter on its own. Pair it with :attr:`consistency_checked`.
        """
        return not np.isfinite(self.q_pvalue) or self.q_pvalue >= DEFAULT_ALPHA

    def ci(self, level: float = 0.95) -> tuple[float, float]:
        """Normal-approximation confidence interval on log2(PTR)."""
        if not np.isfinite(self.se):
            return (np.nan, np.nan)
        z = stats.norm.ppf(0.5 + level / 2.0)
        return (self.log2_ptr - z * self.se, self.log2_ptr + z * self.se)


def _col_min(df: pd.DataFrame, col: str) -> float:
    """Minimum of ``col``, or NaN when the column or the rows are absent."""
    if col not in df.columns or not len(df):
        return float("nan")
    v = pd.to_numeric(df[col], errors="coerce")
    return float(v.min()) if v.notna().any() else float("nan")


def _count_negative(df: pd.DataFrame, col: str) -> int:
    """How many rows have ``col < 0``; 0 when the column is absent."""
    if col not in df.columns or not len(df):
        return 0
    return int((pd.to_numeric(df[col], errors="coerce") < 0).sum())


def fuse(
    estimates: dict[str, float],
    errors: dict[str, float],
    n_anchors: dict[str, int] | None = None,
    alpha: float = DEFAULT_ALPHA,
    min_anchors: int = MIN_ANCHORS_PER_ENZYME,
    random_effects: str = "auto",
) -> FusionResult:
    """Fuse per-enzyme log2(PTR) estimates.

    ``random_effects`` is ``"auto"`` (escalate when Q rejects), ``"always"`` or
    ``"never"``.

    A single surviving enzyme is fused trivially and reported with
    ``n_enzymes = 1``: there is nothing wrong with the estimate, but the
    cross-enzyme QC that motivates this whole design is simply unavailable, and
    the ``note`` says so.
    """
    excluded: dict[str, str] = {}
    names: list[str] = []
    for name in sorted(estimates):
        x, s = estimates[name], errors.get(name, np.nan)
        if not np.isfinite(x):
            excluded[name] = "no estimate"
        elif not np.isfinite(s) or s <= 0:
            excluded[name] = "no usable standard error"
        elif n_anchors is not None and n_anchors.get(name, 0) < min_anchors:
            excluded[name] = f"only {n_anchors.get(name, 0)} anchors (< {min_anchors})"
        else:
            names.append(name)

    if not names:
        return FusionResult(np.nan, np.nan, 0, excluded=excluded, ok=False, note="no enzyme produced a usable estimate")

    x = np.array([estimates[n] for n in names], dtype=float)
    s = np.array([errors[n] for n in names], dtype=float)
    w = 1.0 / s**2

    fixed = float(np.sum(w * x) / np.sum(w))
    q = float(np.sum(w * (x - fixed) ** 2))
    k = len(names)
    dof = k - 1
    if dof > 0:
        pval = float(stats.chi2.sf(q, dof))
        i2 = float(max(0.0, (q - dof) / q)) if q > 0 else 0.0
        # DerSimonian-Laird moment estimator of the between-enzyme variance.
        denom = np.sum(w) - np.sum(w**2) / np.sum(w)
        tau2 = float(max(0.0, (q - dof) / denom)) if denom > 0 else 0.0
    else:
        pval, i2, tau2 = np.nan, np.nan, 0.0

    use_re = random_effects == "always" or (random_effects == "auto" and np.isfinite(pval) and pval < alpha)
    if use_re and tau2 > 0:
        w_eff = 1.0 / (s**2 + tau2)
        est = float(np.sum(w_eff * x) / np.sum(w_eff))
        se = float(np.sqrt(1.0 / np.sum(w_eff)))
        model = "random"
    else:
        est = fixed
        se = float(np.sqrt(1.0 / np.sum(w)))
        model = "fixed"

    residuals = {n: float((xi - est) / si) for n, xi, si in zip(names, x, s)}
    note = ""
    if k == 1:
        note = "single enzyme: no cross-enzyme consistency check available"
    elif np.isfinite(pval) and pval < alpha:
        worst = max(residuals, key=lambda n: abs(residuals[n]))
        note = f"enzymes disagree (Q={q:.1f}, p={pval:.3g}); largest deviation from {worst}"

    return FusionResult(
        log2_ptr=est,
        se=se,
        n_enzymes=k,
        q=q,
        q_pvalue=pval,
        i2=i2,
        tau2=tau2,
        model=model,
        excluded=excluded,
        used=tuple(names),
        residuals=residuals,
        ok=True,
        note=note,
    )


def fuse_table(
    per_enzyme: pd.DataFrame,
    alpha: float = DEFAULT_ALPHA,
    min_anchors: int = MIN_ANCHORS_PER_ENZYME,
    random_effects: str = "auto",
) -> pd.DataFrame:
    """Fuse a per-enzyme fit table into one row per (sample, genome).

    Expects the output of :func:`sk2bgrow.fit.fit_windows`. Rows with
    ``ok == False`` are excluded with their ``note`` recorded, so the reason a
    given enzyme dropped out survives into the report rather than vanishing.
    """
    rows = []
    for (sample, genome_id, genome), grp in per_enzyme.groupby(["sample", "genome_id", "genome"], sort=True):
        est, err, na, pre_excluded = {}, {}, {}, {}
        for _, r in grp.iterrows():
            enzyme = str(r["enzyme"])
            if not bool(r.get("ok", True)):
                pre_excluded[enzyme] = str(r.get("note") or "fit failed")
                continue
            est[enzyme] = float(r["log2_ptr"])
            err[enzyme] = float(r["se"])
            na[enzyme] = int(r.get("n_anchors", 0))
        res = fuse(est, err, na, alpha=alpha, min_anchors=min_anchors, random_effects=random_effects)
        res.excluded.update(pre_excluded)
        n_attempted = int(len(grp))

        used = grp[grp["enzyme"].isin(res.used)]
        lo, hi = res.ci()
        rows.append(
            {
                "sample": sample,
                "genome_id": int(genome_id),
                "genome": genome,
                "log2_ptr": res.log2_ptr,
                "ptr": res.ptr,
                "se": res.se,
                "ci_low": lo,
                "ci_high": hi,
                "n_enzymes": res.n_enzymes,
                "n_enzymes_attempted": n_attempted,
                # An enzyme that failed to fit is itself evidence about the
                # sample — a flat profile, a digestion failure, a mis-assembled
                # region. Without this ratio such enzymes vanish into
                # ``excluded`` and the survivors agree with each other, so a
                # sample where a quarter of the panel saw nothing reads as clean.
                "enzyme_fit_rate": res.n_enzymes / n_attempted if n_attempted else np.nan,
                "enzymes_used": ",".join(res.used) if res.used else "-",
                "enzyme_consistency": res.q_pvalue,
                "enzyme_q": res.q,
                "enzyme_i2": res.i2,
                "tau2": res.tau2,
                "fusion_model": res.model,
                "ok": res.ok,
                "consistent": res.consistent,
                # ``consistent`` is True both when the enzymes agreed and when
                # there was only one of them, so a filter that reads it alone
                # counts unchecked estimates as having passed. This column
                # separates the two cases; require both to mean "checked and
                # agreed".
                "consistency_checked": res.consistency_checked,
                # Per-enzyme fit quality does not otherwise survive fusion.
                # r2 < 0 means the V-fit is worse than a horizontal line, which
                # on the Zheng grid is 57.5% of accepted fits at 0.5x — the
                # ``ok`` flag does not exclude them, so carry the evidence
                # forward instead of dropping it here.
                "min_r2": _col_min(used, "r2"),
                "n_enzymes_negative_r2": _count_negative(used, "r2"),
                "n_anchors": int(used["n_anchors"].sum()) if len(used) else 0,
                "n_windows": int(used["n_windows_used"].sum()) if len(used) else 0,
                "ori": float(used["ori"].median()) if len(used) and used["ori"].notna().any() else np.nan,
                "ori_confidence": float(used["ori_confidence"].median()) if len(used) else np.nan,
                "coverage": float(used["mean_rate"].mean()) if len(used) else np.nan,
                "dispersion": float(used["mean_dispersion"].mean()) if len(used) else np.nan,
                "fraction": float(used["mean_detected_fraction"].mean()) if len(used) else np.nan,
                "method": ",".join(sorted(set(used["method"]))) if len(used) else "none",
                "excluded": ";".join(f"{k}:{v}" for k, v in sorted(res.excluded.items())) or "-",
                "note": res.note,
            }
        )
    return pd.DataFrame(rows)

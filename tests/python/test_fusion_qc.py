"""Executable evidence for what the fused QC does and does not catch.

The C5 run (PRJNA974210, 522 MAGs) reported recall 522/522 = 1.00 in every
sample and zero "suspicious" estimates. Neither number is a property of the
data: the fused output carries no coverage gate, and the flags a consumer would
naturally filter on do not mean what their names suggest. The tests below pin
the parts that are now fixed and mark the parts that are not, so the gap is
visible in CI rather than rediscovered on a cluster.

The ``xfail(strict=True)`` tests assert the behaviour we want. They fail today.
When one is fixed it will XPASS, which strict mode reports as an error — that is
the signal to delete the marker.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from sk2bgrow.fusion import fuse, fuse_table


def enz(
    enzyme: str,
    log2_ptr: float,
    se: float,
    *,
    n_anchors: int = 2000,
    n_windows: int = 29,
    rate: float = 1.5,
    r2: float = 0.7,
    ok: bool = True,
    note: str = "",
) -> dict:
    """One row of the per-enzyme fit table that :func:`fuse_table` consumes."""
    return {
        "sample": "S1",
        "genome_id": 0,
        "genome": "MAG1",
        "enzyme": enzyme,
        "log2_ptr": log2_ptr,
        "ptr": 2.0**log2_ptr if np.isfinite(log2_ptr) else np.nan,
        "se": se,
        "method": "v_shape",
        "n_windows": n_windows,
        "n_windows_used": n_windows,
        "ori": 1_000_000.0,
        "ori_confidence": 0.9,
        "r2": r2,
        "reduced_chi2": 1.0,
        "ok": ok,
        "note": note,
        "mean_rate": rate,
        "mean_detected_fraction": 0.5,
        "mean_dispersion": 0.5,
        "n_anchors": n_anchors,
    }


# --------------------------------------------------------------------------
# Fixed: the fused table now carries the evidence it used to drop.
# --------------------------------------------------------------------------


def test_fused_table_carries_ok():
    """``ok`` survives fusion.

    It did not, so a consumer writing the codebase's own idiom
    ``row.get("ok", True)`` — the pattern :func:`fuse_table` itself uses on the
    per-enzyme table — silently got ``True`` for every genome. That is the same
    defect that once printed a 100% QC pass rate for Pilea in the one table
    about the QC being wrong.
    """
    out = fuse_table(pd.DataFrame([enz("BcgI", 0.4, 0.2), enz("AlfI", 0.45, 0.25)]))
    assert "ok" in out.columns
    assert bool(out["ok"].iloc[0]) is True

    empty = fuse_table(pd.DataFrame([enz("BcgI", np.nan, np.nan, ok=False, note="no gradient")]))
    assert bool(empty["ok"].iloc[0]) is False


def test_single_enzyme_is_consistent_but_not_checked():
    """One enzyme means zero degrees of freedom, so Q never ran.

    ``consistent`` is True — an estimate is not wrong for being alone — but
    ``consistency_checked`` now says the check did not happen, so a QC filter
    can require both.
    """
    out = fuse_table(pd.DataFrame([enz("BcgI", 0.4, 0.2)]))
    assert out["n_enzymes"].iloc[0] == 1
    assert not np.isfinite(out["enzyme_consistency"].iloc[0])
    assert bool(out["consistent"].iloc[0]) is True
    assert bool(out["consistency_checked"].iloc[0]) is False

    two = fuse_table(pd.DataFrame([enz("BcgI", 0.4, 0.2), enz("AlfI", 0.42, 0.22)]))
    assert bool(two["consistency_checked"].iloc[0]) is True


def test_negative_r2_is_carried_forward():
    """A V-fit with r2 < 0 is worse than a horizontal line.

    On the Zheng grid 15.9% of accepted per-enzyme fits are in that state, and
    57.5% at 0.5x. Fusion used to drop ``r2`` entirely, so nothing downstream
    could see it.
    """
    out = fuse_table(
        pd.DataFrame([enz("BcgI", 0.4, 0.2, r2=-0.03), enz("AlfI", 0.45, 0.25, r2=0.6)])
    )
    assert out["min_r2"].iloc[0] == pytest.approx(-0.03)
    assert int(out["n_enzymes_negative_r2"].iloc[0]) == 1


# --------------------------------------------------------------------------
# Not fixed: each of these changes reported numbers, so they wait for the
# C5 rescore rather than landing mid-run.
# --------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="no coverage gate exists anywhere in fit or fuse")
def test_trace_coverage_estimate_is_withheld():
    """A genome at 0.02 counts per anchor should not receive a PTR estimate.

    This is the mechanism behind C5's recall = 1.00 and suspicious = 0: the
    pipeline has no equivalent of Pilea's coverage gate, so it reports for every
    genome in the index whether or not the sample contains it.
    """
    out = fuse_table(pd.DataFrame([enz("BcgI", 0.02, 0.30, n_windows=6, rate=0.02)]))
    assert bool(out["ok"].iloc[0]) is False


@pytest.mark.xfail(strict=True, reason="coverage averages surviving enzymes only, unweighted")
def test_reported_coverage_is_not_survivorship_biased():
    """``coverage`` must describe the genome, not the enzymes that succeeded.

    Enzymes drop out *because* coverage was too low, so averaging the survivors
    reports a number systematically above the truth — 6.5x in this case. Any
    "suspicious" detector thresholding on this column will not fire.
    """
    rows = [enz("BcgI", 0.4, 0.2, rate=1.6), enz("AlfI", 0.45, 0.25, rate=1.5)]
    rows += [
        enz(f"E{i}", np.nan, np.nan, n_anchors=300, n_windows=3, rate=0.05, ok=False,
            note="fewer than 5 usable windows")
        for i in range(14)
    ]
    truth = float(np.mean([r["mean_rate"] for r in rows]))
    out = fuse_table(pd.DataFrame(rows))
    assert out["coverage"].iloc[0] == pytest.approx(truth, rel=0.25)


@pytest.mark.xfail(strict=True, reason="fuse ignores r2; fit's ok flag is only a sign check")
def test_worse_than_constant_fit_is_excluded_from_fusion():
    """An enzyme whose fit is worse than a constant should not be fused in.

    ``fit.py`` sets ``ok = log2_ptr >= 0``, which on 1,360 real fits rejected
    nothing — all 87 rejections came from earlier return paths, and zero had a
    negative estimate. So the flag is a sign check, not a quality check, and
    ``fuse`` never looks at ``r2`` at all.
    """
    res = fuse(
        estimates={"BcgI": 0.40, "AlfI": 0.45},
        errors={"BcgI": 0.20, "AlfI": 0.25},
        n_anchors={"BcgI": 2000, "AlfI": 1800},
        r2={"BcgI": -0.03, "AlfI": 0.60},  # type: ignore[call-arg]
    )
    assert res.n_enzymes == 1
    assert "BcgI" in res.excluded

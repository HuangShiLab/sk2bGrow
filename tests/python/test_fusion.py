"""Cross-enzyme fusion — the project's core new step."""

import numpy as np
import pandas as pd
import pytest

from sk2bgrow import fusion


def enzyme_set(values, se=0.05, n=500):
    est = {f"E{i:02d}": v for i, v in enumerate(values)}
    err = {k: se for k in est}
    cnt = {k: n for k in est}
    return est, err, cnt


def test_fusion_reduces_variance_by_root_n():
    """16 independent enzymes at se=0.05 must fuse to se=0.0125."""
    est, err, cnt = enzyme_set([1.0] * 16)
    r = fusion.fuse(est, err, cnt)
    assert r.n_enzymes == 16
    assert r.se == pytest.approx(0.05 / 4.0, rel=1e-6)
    assert r.log2_ptr == pytest.approx(1.0, abs=1e-9)
    assert r.model == "fixed"


def test_precise_enzymes_get_more_weight():
    est = {"A": 1.0, "B": 2.0}
    err = {"A": 0.1, "B": 0.2}  # A carries 4x the information of B
    # Pinned to the fixed-effect model so this checks the weighting formula
    # itself; these two estimates are heterogeneous enough that "auto" would
    # (correctly) escalate to random effects and flatten the weights.
    r = fusion.fuse(est, err, {"A": 500, "B": 500}, random_effects="never")
    assert r.log2_ptr == pytest.approx((1.0 / 0.01 + 2.0 / 0.04) / (1 / 0.01 + 1 / 0.04))
    assert r.log2_ptr == pytest.approx(1.2, abs=1e-9)


def test_random_effects_flatten_weights_when_enzymes_disagree():
    est = {"A": 1.0, "B": 2.0}
    err = {"A": 0.1, "B": 0.2}
    fixed = fusion.fuse(est, err, {"A": 500, "B": 500}, random_effects="never")
    auto = fusion.fuse(est, err, {"A": 500, "B": 500})
    # With a large between-enzyme variance the precision difference matters
    # less, so the estimate moves toward the unweighted mean of 1.5.
    assert auto.model == "random"
    assert fixed.log2_ptr < auto.log2_ptr <= 1.5


def test_consistent_enzymes_pass_the_q_test():
    rng = np.random.default_rng(0)
    est, err, cnt = enzyme_set(1.0 + rng.normal(0, 0.05, 16))
    r = fusion.fuse(est, err, cnt)
    assert r.consistent
    assert r.q_pvalue > 0.05
    assert r.i2 < 0.5


def test_one_rogue_enzyme_is_detected_and_named():
    values = [1.0] * 16
    values[3] = 2.5
    est, err, cnt = enzyme_set(values)
    r = fusion.fuse(est, err, cnt)
    assert not r.consistent
    assert r.q_pvalue < 1e-6
    assert r.i2 > 0.9
    assert "E03" in r.note
    assert abs(r.residuals["E03"]) == max(abs(v) for v in r.residuals.values())


def test_heterogeneity_widens_the_interval():
    """A tight interval around a value the enzymes disagree about would be a lie."""
    est_c, err, cnt = enzyme_set([1.0] * 16)
    values = [1.0] * 16
    values[3] = 2.5
    est_h, _, _ = enzyme_set(values)
    consistent = fusion.fuse(est_c, err, cnt)
    heterogeneous = fusion.fuse(est_h, err, cnt)
    assert heterogeneous.model == "random"
    assert heterogeneous.se > 5 * consistent.se
    assert heterogeneous.tau2 > 0


def test_random_effects_can_be_disabled():
    values = [1.0] * 16
    values[3] = 2.5
    est, err, cnt = enzyme_set(values)
    r = fusion.fuse(est, err, cnt, random_effects="never")
    assert r.model == "fixed"
    assert r.se == pytest.approx(0.0125, rel=1e-6)
    assert not r.consistent, "the Q test must still fire even without escalation"


def test_thin_enzymes_are_excluded_with_a_reason():
    est, err, cnt = enzyme_set([1.0] * 4)
    cnt["E01"] = 5
    r = fusion.fuse(est, err, cnt, min_anchors=30)
    assert "E01" not in r.used
    assert "5 anchors" in r.excluded["E01"]
    assert r.n_enzymes == 3


def test_unusable_estimates_are_excluded():
    est = {"A": 1.0, "B": np.nan, "C": 1.0}
    err = {"A": 0.1, "B": 0.1, "C": np.nan}
    r = fusion.fuse(est, err, {"A": 500, "B": 500, "C": 500})
    assert r.used == ("A",)
    assert r.excluded == {"B": "no estimate", "C": "no usable standard error"}
    assert "single enzyme" in r.note


def test_no_usable_enzyme_fails_loudly():
    r = fusion.fuse({"A": np.nan}, {"A": np.nan}, {"A": 0})
    assert not r.ok
    assert np.isnan(r.log2_ptr)


def test_single_enzyme_has_no_consistency_check():
    r = fusion.fuse({"A": 1.0}, {"A": 0.1}, {"A": 500})
    assert r.n_enzymes == 1
    assert np.isnan(r.q_pvalue)
    assert r.consistent, "an untestable hypothesis is not a failed one"
    assert "no cross-enzyme consistency check" in r.note


def test_confidence_interval_covers_the_estimate():
    r = fusion.fuse(*enzyme_set([1.0] * 16))
    lo, hi = r.ci(0.95)
    assert lo < r.log2_ptr < hi
    assert hi - lo == pytest.approx(2 * 1.959964 * r.se, rel=1e-4)


def test_fuse_table_groups_by_sample_and_genome():
    rows = []
    for sample in ["S1", "S2"]:
        for i, enzyme in enumerate(["BcgI", "AlfI", "CjeI"]):
            rows.append(
                {
                    "sample": sample, "genome_id": 0, "genome": "ecoli", "enzyme": enzyme,
                    "log2_ptr": 1.0 + 0.01 * i, "se": 0.05, "n_anchors": 500,
                    "n_windows_used": 12, "ori": 3_923_883, "ori_confidence": 0.9,
                    "mean_rate": 8.0, "mean_dispersion": 1.1, "mean_detected_fraction": 0.98,
                    "method": "v_shape", "ok": True, "note": "",
                }
            )
    out = fusion.fuse_table(pd.DataFrame(rows))
    assert len(out) == 2
    assert set(out["sample"]) == {"S1", "S2"}
    assert (out["n_enzymes"] == 3).all()
    assert (out["consistent"]).all()
    assert out["enzymes_used"].iloc[0] == "AlfI,BcgI,CjeI"


def test_fuse_table_carries_failed_fits_into_excluded():
    rows = [
        {"sample": "S1", "genome_id": 0, "genome": "g", "enzyme": "BcgI", "log2_ptr": 1.0, "se": 0.05,
         "n_anchors": 500, "n_windows_used": 10, "ori": 1, "ori_confidence": 1.0, "mean_rate": 5,
         "mean_dispersion": 1, "mean_detected_fraction": 1, "method": "v_shape", "ok": True, "note": ""},
        {"sample": "S1", "genome_id": 0, "genome": "g", "enzyme": "PpiI", "log2_ptr": np.nan, "se": np.nan,
         "n_anchors": 12, "n_windows_used": 1, "ori": np.nan, "ori_confidence": np.nan, "mean_rate": np.nan,
         "mean_dispersion": np.nan, "mean_detected_fraction": np.nan, "method": "none", "ok": False,
         "note": "only 1 windows after outlier trimming"},
    ]
    out = fusion.fuse_table(pd.DataFrame(rows))
    assert out["n_enzymes"].iloc[0] == 1
    assert "PpiI" in out["excluded"].iloc[0]
    assert "only 1 windows" in out["excluded"].iloc[0]

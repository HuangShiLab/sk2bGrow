"""Multi-sample delta-PTR."""

import numpy as np
import pandas as pd
import pytest

from sk2bgrow import dynamics


def outputs(rows):
    return pd.DataFrame(
        [
            {"sample": s, "genome": g, "log2(PTR)": v, "se": se, "pass_qc": qc}
            for s, g, v, se, qc in rows
        ]
    )


def test_delta_against_a_named_baseline():
    df = outputs([("T0", "ecoli", 1.0, 0.05, True), ("T1", "ecoli", 2.0, 0.05, True)])
    d = dynamics.delta_ptr(df, baseline="T0")
    assert len(d) == 1
    r = d.iloc[0]
    assert r["sample"] == "T1"
    assert r["delta_log2_ptr"] == pytest.approx(1.0)
    assert r["fold_change"] == pytest.approx(2.0)
    # Independent estimates: variances add.
    assert r["se"] == pytest.approx(np.sqrt(0.05**2 + 0.05**2))
    assert r["p"] < 1e-30


def test_no_baseline_uses_the_per_genome_mean():
    df = outputs([("A", "g", 1.0, 0.1, True), ("B", "g", 2.0, 0.1, True), ("C", "g", 3.0, 0.1, True)])
    d = dynamics.delta_ptr(df)
    assert set(d["sample"]) == {"A", "B", "C"}
    assert d.loc[d["sample"] == "B", "delta_log2_ptr"].iloc[0] == pytest.approx(0.0)
    assert d["baseline_log2_ptr"].nunique() == 1


def test_qc_failures_are_excluded_by_default():
    df = outputs([("T0", "g", 1.0, 0.05, True), ("T1", "g", 9.0, 0.05, False)])
    assert dynamics.delta_ptr(df, baseline="T0").empty
    kept = dynamics.delta_ptr(df, baseline="T0", qc_only=False)
    assert len(kept) == 1


def test_genome_absent_from_the_baseline_is_skipped():
    df = outputs([("T0", "a", 1.0, 0.05, True), ("T1", "a", 2.0, 0.05, True), ("T1", "b", 3.0, 0.05, True)])
    d = dynamics.delta_ptr(df, baseline="T0")
    assert set(d["genome"]) == {"a"}, "genome b has no baseline and cannot be differenced"


def test_group_baseline_from_metadata():
    df = outputs([("c1", "g", 1.0, 0.05, True), ("c2", "g", 1.1, 0.05, True), ("t1", "g", 2.0, 0.05, True)])
    meta = pd.DataFrame({"sample": ["c1", "c2", "t1"], "group": ["ctrl", "ctrl", "treat"]})
    d = dynamics.delta_ptr(df, baseline="ctrl", metadata=meta)
    assert list(d["sample"]) == ["t1"]
    assert d["delta_log2_ptr"].iloc[0] == pytest.approx(2.0 - 1.05)


def test_bh_correction_is_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.5, np.nan])
    q = dynamics._bh(p)
    assert np.isnan(q[3])
    assert np.all(q[:3] >= p[:3])
    assert np.all(np.diff(q[:3]) >= -1e-12)
    assert np.all(q[:3] <= 1.0)


def test_trend_test_detects_a_rising_series():
    rows = []
    for t, v in enumerate([1.0, 1.5, 2.0, 2.5, 3.0]):
        rows.append({"genome": "g", "sample": f"T{t}", "timepoint": t, "log2_ptr": v, "se": 0.05})
    tr = dynamics.trend_test(pd.DataFrame(rows))
    assert tr["slope"].iloc[0] == pytest.approx(0.5, abs=0.02)
    assert tr["p"].iloc[0] < 0.01


def test_trend_test_needs_three_timepoints():
    rows = [{"genome": "g", "sample": "T0", "timepoint": 0, "log2_ptr": 1.0, "se": 0.1},
            {"genome": "g", "sample": "T1", "timepoint": 1, "log2_ptr": 2.0, "se": 0.1}]
    tr = dynamics.trend_test(pd.DataFrame(rows))
    assert np.isnan(tr["slope"].iloc[0])
    assert "3 timepoints" in tr["note"].iloc[0]


def test_anchor_matrix_keys_are_stable_across_samples(tmp_path):
    from conftest import make_counts
    from sk2bgrow import io as sk_io

    a = make_counts(n_per_enzyme=60, enzymes=("BcgI",), sample="S1", seed=0)
    b = a.copy()
    b["sample"] = "S2"
    b["count"] = b["count"] * 2
    paths = []
    for name, frame in [("S1", a), ("S2", b)]:
        p = tmp_path / f"{name}.counts.tsv"
        frame[sk_io.COUNT_COLUMNS].to_csv(p, sep="\t", index=False)
        paths.append(p)
    m = dynamics.anchor_matrix(paths)
    # The deterministic sketch is what makes this matrix well defined: the same
    # loci are observed in both samples, so there are no missing cells.
    assert list(m.columns) == ["S1", "S2"]
    assert m.notna().all().all()
    assert (m["S2"] >= m["S1"]).all()


def test_ptr_matrix_pivots_genomes_by_sample():
    df = outputs([("S1", "a", 1.0, 0.05, True), ("S2", "a", 2.0, 0.05, True), ("S1", "b", 0.5, 0.05, True)])
    m = dynamics.ptr_matrix(df)
    assert m.loc["a", "S2"] == 2.0
    assert np.isnan(m.loc["b", "S2"])
